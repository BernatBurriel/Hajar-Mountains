# To run this script you will need several files: 
	# genotype file in .bim, .bam .fam  format
	# decimal coordinates of the samples 
	# coordinates of the outer margin of the distribution
	# shapefile containing a ddg grid
# base 
import math
import numpy as np
import pkg_resources
from sklearn.impute import SimpleImputer
from pandas_plink import read_plink
from sklearn.metrics.pairwise import haversine_distances
from scipy.spatial.distance import pdist, squareform
import statsmodels.api as sm
import pickle

# export needed module
import fiona

# viz
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib import gridspec

# feems
from feems.utils import prepare_graph_inputs
from feems import SpatialGraph, Viz, Objective, query_node_attributes
from feems.cross_validation import run_cv, comp_mats



# Extra code to add export function to Viz

def to_shapefile(self, filename):
    """ Exports a grid with attributes to a line shapefile"""
    crs = self.projection.proj4_params
    schema = {'geometry': 'LineString',
              'properties': {
                'feems': 'float'}}
    with fiona.open(
        filename, "w", 
        crs=crs, 
        driver="ESRI Shapefile", 
        schema=schema) as stream:
        
        e_done = []
        for i in range(len(self.weights)):
            #Checks duplicates by sorting the indices and testing against e_done
            e = list(np.sort(np.column_stack(self.idx)[i]))
            if e not in e_done:
                l = {"type": "LineString",
                     "coordinates": self.grid[e].tolist()}
                p = {'feems': self.norm_log_weights[i] }
                
                stream.write({'geometry':l, 'properties': p})
                e_done.append(e)


Viz.to_shapefile = to_shapefile


## END of EXTRA code ##



# GIVE PATH TO THE DATA
data_path2 = "/home/pristurus/Desktop/Bernat/feems/DATA/TEST_FEEMS/"

# read the genotype data and mean impute missing data
(bim, fam, G) = read_plink("{}/PATH_TO_FILE".format(data_path2))
imp = SimpleImputer(missing_values=np.nan, strategy="mean")
genotypes = imp.fit_transform((np.array(G)).T)

print("n_samples={}, n_snps={}".format(genotypes.shape[0], genotypes.shape[1]))

# setup graph 
coord = np.loadtxt("{}/PATH_TO_FILE.coords".format(data_path2))  # sample coordinates
outer = np.loadtxt("{}/outer_coordinates.outer".format(data_path2))  # outer coordinates
grid_path = "{}/shapefile_to_grid_10km.shp".format(data_path2)  # path to discrete global grid

# graph input files
outer, edges, grid, _ = prepare_graph_inputs(coord=coord, 
                                             ggrid=grid_path,
                                             translated=False, 
                                             buffer=0,
                                             outer=outer)

# We then setup the SpatialGraph object which is the core workhorse of feems. SpatialGraph specifies the graph, allele frequency data, and runs 
# the optimizers to fit the edge weights of the graph:

sp_graph = SpatialGraph(genotypes, coord, grid, edges, scale_snps=True)

#  This might take a few minutes to construct at first b/c it initializing a number of graph matrices that are slow to build. 
# First, before any fitting we’ll visualize the graph and samples. Let’s set up the projection we’ll be using for this dataset:

projection = ccrs.EquidistantConic(central_longitude= 57.27097 , central_latitude=23.44036 ) # fixed into the center of the Hajars buffer30km

# plot several lambdas
# fit
sp_graph.fit(100.0)

# Prepare plot
fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=20, 
        obs_node_size=7.5, sample_pt_color="black", 
        cbar_font_size=10)
        
# export to shapefile
v.to_shapefile("FILENAME100.shp")

# lambda 10
sp_graph.fit(10.0)

# Prepare plot
fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=20, 
        obs_node_size=7.5, sample_pt_color="black", 
        cbar_font_size=10)
        
# export to shapefile
v.to_shapefile("FILENAME10.shp")

#lambda 1
sp_graph.fit(1.0)

# Prepare plot
fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=20, 
        obs_node_size=7.5, sample_pt_color="black", 
        cbar_font_size=10)
        
# export to shapefile
v.to_shapefile("FILENAME1.shp")

#lamda 0.1
sp_graph.fit(0.1)

# Prepare plot
fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=20, 
        obs_node_size=7.5, sample_pt_color="black", 
        cbar_font_size=10)
        
# export to shapefile
v.to_shapefile("FILENAME01.shp")

#lambda 0.01
sp_graph.fit(0.01)

# Prepare plot
fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=20, 
        obs_node_size=7.5, sample_pt_color="black", 
        cbar_font_size=10)
        
# export to shapefile
v.to_shapefile("FILENAME0.01.shp")

# lambda 0.001
sp_graph.fit(0.001)

# Prepare plot
fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=20, 
        obs_node_size=7.5, sample_pt_color="black", 
        cbar_font_size=10)
        
# export to shapefile
v.to_shapefile("FILENAME0.001.shp")

#### Run cross-validation
# Next we perform leave-one-out cross-valiation over a grid of $\lambda$ values. In our CV we hold out indivudal observed nodes on the graph, predict allele frequencies 
# at the held-out node under our fitted model, and use the $\ell_2$ distance between the fitted and predicted allele frequencies as our CV metric to select models:

# define grids
# reverse the order of lambdas and alphas for warmstart
lamb_grid = np.geomspace(1e-6, 1e2, 20)[::-1]

# run cross-validation
cv_err = run_cv(sp_graph, lamb_grid, n_folds=sp_graph.n_observed_nodes, factr=1e10)

# average over folds
mean_cv_err = np.mean(cv_err, axis=0)

# argmin of cv error
lamb_cv = float(lamb_grid[np.argmin(mean_cv_err)])
lamb_cv

# re-fit
sp_graph.fit(lamb_cv)

# Prepare plot
fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=20, 
        obs_node_size=7.5, sample_pt_color="black", 
        cbar_font_size=10)
        
# export to shapefile
v.to_shapefile("FILENAME.bestCV.shp")
