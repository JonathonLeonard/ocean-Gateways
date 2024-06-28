#!/usr/bin/python3

import astropy.units
from fil_finder import FilFinder2D
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import os
import pygmt
import shutil
from skimage.morphology import skeletonize, flood_fill
import sys

def generate_gateway(grid, ocean_or_land='ocean'):

	outname = f'skeleton_{ocean_or_land}_{grid}'
	if os.path.exists(outname):
		os.remove(outname)

	shutil.copy(f'{grid}', os.path.join(os.getcwd(), outname))
	grid = outname

	# Translate grid to 0-360 lon and the pixel registration
	# pygmt.grdsample(grid=grid, region='d', registration='p', outgrid=grid)

	# Open the netCDF file
	dataset = nc.Dataset(grid, mode='r+')

	# Extract the variable as a NumPy array
	grid_data = np.array(dataset.variables['landsea_smoothed_nosea_noisland_nodiagonal'])

	grid_data = grid_data.astype(int) # Make every value an integer

	if ocean_or_land == 'land':
		# Create land and water masks
		# print('nothing')
		pass
	elif ocean_or_land == 'ocean':
		# Use logical_not to reverse the array
		grid_data = np.logical_not(grid_data).astype(int)

	# Get a basic skelton using skimage
	skel = skeletonize(grid_data)

	# Create a new variable in the dataset to store the 'skel' array
	var = dataset.createVariable('skel', 'i4', dataset.variables['z'].dimensions)

	# Save the 'skel' array to the new variable
	var[:] = skel

	# Process and prune the skeleton with Fil_Finder
	# Note that this is adapted from the online doc, it's possible not all these steps are needed
	fil = FilFinder2D(skel, mask=skel)
	fil.create_mask(border_masking=True, verbose=False, use_existing_mask=True)
	fil.medskel(verbose=False)
	# Change branch thresh below to change the amount of pruning
	fil.analyze_skeletons(skel_thresh=3 * astropy.units.pix, branch_thresh=100* astropy.units.pix, prune_criteria='length')

	# Create a new variable in the dataset to store the 'skel' array
	var = dataset.createVariable('pruned_skel', 'i4', dataset.variables['z'].dimensions)

	# Save the 'skel' array to the new variable
	var[:] = fil.skeleton

	# Trim the ends of the skeleton so they don't mess with coastlines and other skeletons etc
	skeleton = fil.skeleton
	loose_ends = []
	iteration = 0
	while iteration <= 20: # Trims the ends by ~20 cells- can adjust this as needed
		for i in range(len(skeleton) - 1):
			for j in range(len(skeleton[i]) - 1):
				cell = skeleton[i][j]
				if not cell: # Only pay attention to the skeleton cells (the cells that are True)
					continue

				# Add the surrounding cells, their value and the indicies to a dictionary. Also ignore if the neighbouring cells is over the edge
				connections = {
					'above': (skeleton[i-1][j], (i-1, j)) if skeleton[i-1][j] else (False, None),
					'below': (skeleton[i+1][j], (i+1, j)) if skeleton[i+1][j] else (False, None),
					'left': (skeleton[i][j-1], (i, j-1)) if skeleton[i][j-1] else (False, None),
					'right': (skeleton[i][j+1], (i, j+1)) if skeleton[i][j+1] else (False, None),
					'upleft': (skeleton[i-1][j-1], (i-1, j-1)) if skeleton[i-1][j-1] else (False, None),
					'upright': (skeleton[i-1][j+1], (i-1, j+1)) if skeleton[i-1][j+1] else (False, None),
					'lowleft': (skeleton[i+1][j-1], (i+1, j-1)) if skeleton[i+1][j-1] else (False, None),
					'lowright': (skeleton[i+1][j+1], (i+1, j+1)) if skeleton[i+1][j+1] else (False, None)
				}

				# Trim if there's only 1 connection (or 0, i.e. it's a standalone skeleton cell) i.e. it's at the end of the skeleton
				if sum(value for (value, indices) in connections.values()) <= 1:
					loose_ends.append((i,j))

					
		for (i,j) in loose_ends:
			skeleton[i][j] = False
		iteration += 1

	# Create a new variable in the dataset to store the 'skel' array
	var = dataset.createVariable('pruned_trimmed_skel', 'i4', dataset.variables['z'].dimensions)

	# Save the 'skel' array to the new variable
	var[:] = skeleton

	# Edit the downscaled skeleton to give diagonal cells an adjacent buddy
	for i in range(len(skeleton) - 1):
		for j in range(len(skeleton[i])):
			
			if j == len(skeleton[i]) - 1:
				continue

			cell = skeleton[i][j]
			if not cell:
				continue

			above = skeleton[i-1][j]
			below = skeleton[i+1][j]
			left = skeleton[i][j-1]
			right = skeleton[i][j+1]
			upleft = skeleton[i-1][j-1]
			upright = skeleton[i-1][j+1]
			lowleft = skeleton[i+1][j-1]
			lowright = skeleton[i+1][j+1]

			
			if (above or below) and (right or left):
				continue

			if upleft and not (left or above):
				skeleton[i-1][j] = True  # modify the 'above' cell

			if lowleft and not (left or below):
				skeleton[i+1][j] = True  # modify the 'below' cell

			if lowright and not (right or below):
				skeleton[i+1][j] = True  # modify the 'below' cell
			
			if upright and not (right or above):
				skeleton[i-1][j] = True  # modify the 'above' cell

	# Create a new variable in the dataset to store the 'skel' array
	var = dataset.createVariable('pruned_trimmed_skel_noDiagonal', 'i4', dataset.variables['z'].dimensions)

	# Save the 'skel' array to the new variable
	var[:] = skeleton

	# Close the dataset
	dataset.close()

	return grid

def makeSmoothMask(grid):
	"""
	Takes a netCDF grid of elevations and produces a smoothed landsea mask
	Use for making skeletons
	"""

	# Open the netCDF file
	dataset = nc.Dataset(grid, mode='r+')

	# Extract the variable as a NumPy array
	grid_data = np.array(dataset.variables['z']).astype(int)

	land = np.where(grid_data > 0)
	water = np.where(grid_data <= 0)

	grid_data[land] = 0 # Isolate all land elevations and make it a high number
	grid_data[water] = 1 # Make all water cells zero

	# Create a new variable in the dataset to store the 'skel' array
	landsea_var = dataset.createVariable('landsea', 'i4', dataset.variables['z'].dimensions)

	# Save the 'skel' array to the new variable
	landsea_var[:] = grid_data

	nosea_smoothed_grid_data = remove_inland_seas(grid_data, minimum_sea_size=500)

	# Create a new variable in the dataset to store the 'skel' array
	var = dataset.createVariable('landsea_smoothed_nosea', 'i4', dataset.variables['z'].dimensions)

	# Save the smoothed array to the new variable
	var[:] = nosea_smoothed_grid_data

	nosea_noisland_smoothed_grid_data = remove_islands(nosea_smoothed_grid_data, minimum_island_size=500)

	# Put the grid back in the form we had originally
	nosea_noisland_smoothed_grid_data[nosea_noisland_smoothed_grid_data == 0] = 1
	nosea_noisland_smoothed_grid_data[nosea_noisland_smoothed_grid_data == 999] = 0

	# plt.matshow(grid_data)
	# plt.matshow(nosea_noisland_smoothed_grid_data)
	# sys.exit()

	# Create a new variable in the dataset to store the 'skel' array
	var = dataset.createVariable('landsea_smoothed_nosea_noisland', 'i4', dataset.variables['z'].dimensions)

	# Save the smoothed array to the new variable
	var[:] = nosea_noisland_smoothed_grid_data

	finalArray = remove_diagonal_connections(nosea_noisland_smoothed_grid_data)

	# Create a new variable in the dataset to store the 'skel' array
	var = dataset.createVariable('landsea_smoothed_nosea_noisland_nodiagonal', 'i4', dataset.variables['z'].dimensions)

	# Save the smoothed array to the new variable
	var[:] = nosea_noisland_smoothed_grid_data

	dataset.close()

	return grid

def remove_inland_seas(grid, minimum_sea_size):
	"""Uses a flood fill algorithm to identify unique bodies of water and removes any below the min waterbody size (set in config)"""
	grid = grid.astype(int) # Make every value an integer

	print(f'Identifying unique bodies of water')

	# Create land and water masks
	land = np.where(grid > 0)
	water = np.where(grid <= 0)

	grid[land] = 999 # Isolate all land elevations and make it a high number
	grid[water] = 0 # Make all water cells zero

	# Initialise variables
	unfilled = np.where(grid == 0) # Everywhere that the flood_fill algorithm hasn't filled yet. (Will be all ocean for now)
	seed_lats = unfilled[0] # A list of all the y-coords that are unfilled
	seed_lons = unfilled[1] # A list of all the x-coords thart are unfilled

	# Run the flood fill algorithm for the first seed lat and seed lon on the list. Everything connected to that (connectiveity=1 doesn't do diagonals) will be a 1
	# print(f"Doing flood fill #1 at {seed_lats[0]},{seed_lons[0]}")
	filled = flood_fill(grid, (seed_lats[0],seed_lons[0]), 1, connectivity=1)

	# Now iterate over unfilled areas until there are no unfilled areas left
	i=2
	while len(np.where(filled == 0)[0]) > 0:
		unfilled = np.where(filled == 0)
		seed_lats = unfilled[0]
		seed_lons = unfilled[1]
		# print(f'Doing floodfill #{i} at {seed_lats[0]},{seed_lons[0]}')

		filled = flood_fill(filled, (seed_lats[0],seed_lons[0]), i, connectivity=1)
		i += 1

	# Check we haven't labelled connected water across the edges of the map as unconnected
	print(f"Floodfill of water bodies done")
	print(f"Now checking connectivity of water across grid edges")
	for r in range(len(filled)):
		# print(f' Floodfill: Left edge is {filled[r][0]}, right edge is {filled[r][-1]}')
		if (filled[r][0] < 999 and filled[r][-1] < 999) and (filled[r][0] != filled[r][-1]):
			print('uh oh, there is connected water across the grid edges')
			print('fixing..')
			filled = flood_fill(filled, (r,-1), filled[r][0], connectivity=1)

	uniques = list(np.unique(filled))
	uniques.remove(999)
	print(f'Detected {len(uniques)} unique bodies of water')

	count = 1
	for sea in uniques:
		print(f'Looking at water body #{count}')
		size = np.count_nonzero(filled == sea)
		print(f'This water body takes up {size} cells')
		if size < minimum_sea_size:
			print(f'Min sea size is {minimum_sea_size}, so turning this to land')
			mask = np.where(filled == sea)
			filled[mask] = 999
		else:
			print('This water body is big enough to stay')
		count += 1

	return filled

def remove_islands(grid,minimum_island_size):
	"""Removes islands smaller than the min island size (see config)"""

	print(f'Identifying unique islands')

	# Create land and water masks
	land = np.where(grid == 999)
	water = np.where(grid < 999)

	grid[water] = 0 # Now that water bodies are filters, we return water to zero

	seed_lats = land[0]
	seed_lons = land[1]
	print(f"Doing flood fill #1 at {seed_lats[0]},{seed_lons[0]}")
	filled = flood_fill(grid, (seed_lats[0],seed_lons[0]), 1, connectivity=1)

	i = 2
	while len(np.where(filled == 999)[0]) > 0:
		unfilled = np.where(filled == 999)
		seed_lats = unfilled[0]
		seed_lons = unfilled[1]
		print(f"Doing flood fill #{i} at {seed_lats[0]},{seed_lons[0]}")
		filled = flood_fill(filled, (seed_lats[0],seed_lons[0]), i, connectivity=1)

		i += 1

	# Check we haven't labelled connected land across the edges of the map as unconnected
	for r in range(len(filled)): # loop over row indicies
		print(f' Floodfill islands: Left edge is {filled[r][0]}, right edge is {filled[r][-1]}')
		if (filled[r][0] != 0 and filled[r][-1] != 0) and (filled[r][0] != filled[r][-1]):
			print('uh oh, there is connected water across the grid edges')
			print('putting in code to fix this but has not been tested yet, so better double check')

			print(f'seed {r,-1}')
			print(f'fill value {filled[r][0]}')

			filled = flood_fill(filled, (r,-1), filled[r][0], connectivity=1)

	# Now convert any islands that are too small to ocean
	uniques = list(np.unique(filled))
	uniques.remove(0) # 0 is water so we don't count this as an island

	print(f'Detected {len(uniques)} unique islands')

	count=1
	for island in uniques:
		print(f'Looking at island #{count}')
		size = np.count_nonzero(filled == island)
		print(f'This island takes up {size} cells')
		if size < minimum_island_size:
			print(f'Min island size is {minimum_island_size}, so turning this to land')
			mask = np.where(filled == island)
			filled[mask] = 0
		else:
			print('This island is big enough to stay')
		count += 1	

	# Convert land back to the 999 value
	filled[np.where(filled > 0)] = 999	

	return filled

def remove_diagonal_connections(array):
	""" 
	Removes diagonal-only connections as this will create a 'crossing' of skeletons
	"""
	for i in range(len(array) - 1):
		for j in range(len(array[i])):
			if j == len(array[i]) - 1:
				continue

			cell = array[i][j]
			if cell == 1: # Only pay attention to land connections (cells that equal 0)
				continue

			above = array[i-1][j]
			below = array[i+1][j]
			left = array[i][j-1]
			right = array[i][j+1]
			upleft = array[i-1][j-1]
			upright = array[i-1][j+1]
			lowleft = array[i+1][j-1]
			lowright = array[i+1][j+1]

			if (upleft == 0 and above == 1 and left == 1) or (upright == 0 and above == 1 and right == 1):
				print(i,j)
				print(f'{upleft} {above} {upright}')
				print(f'{left} {cell} {right}')
				print(f'{lowleft} {below} {lowright}')
				# array[i-1][j] = 0 # For if you want the cell to become a landbridge
				array[i][j] = 1 # For if you want the area to become an ocean gateway

			elif (lowright == 0 and below == 1 and right == 1) or (lowleft == 0 and below == 1 and left == 1):
				print(i,j)
				print(f'{upleft} {above} {upright}')
				print(f'{left} {cell} {right}')
				print(f'{lowleft} {below} {lowright}')
				# array[i+1][j] = 0 # For if you want the cell to become a landbridge
				array[i][j] = 1 # For if you want the area to become an ocean gateway
		return array
