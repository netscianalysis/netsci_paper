# This is designed to be used for analysis with other scripts or jupyter notebooks. 
# It may be converted into a script that can be used frequently from the command line. 

### Updated by Marcus Hock 05/23/2024


import matplotlib.pyplot as plt
import numpy as np
import mdtraj as md
import networkx as nx
import seaborn as sns
from tqdm import tqdm 
import mdtraj as md

color_scheme = ['#777777','#332288', '#44AA99']

apo_color = color_scheme[0]
ATP_color = color_scheme[1]
dATP_color = color_scheme[2]

# Funciton that multiplies or uses a cutoff. 
# Cutoff should be in nm (likely 1.5 nm)
# Note: Topology file should ideally be PDB and the coordinates that are desired to be used 
# for the contact map, unless a separate argument is to be used. 
def correlation_to_distance(correlation_array, topology, cutoff = None, input_contact_dist = None):
    # Define new adjacency matrix to return 

    # Attempt to load it directly as the file type. 
    if type(correlation_array) != np.ndarray:
        try: 
            temp_corr_array = np.load(correlation_array)
            print('Loaded correlation as nd array.')
            
        except:
            print('Unable to load numpy nd array. Trying as text load. ')

            try: 
                temp_corr_array = np.loadtxt(correlation_array)
                print('Loaded correlation as text.')
            except:
                print('Unable to load as text either, exiting.')
                exit()
        correlation_array = temp_corr_array


    log_dist = -np.log((correlation_array))
    log_dist[log_dist == np.inf] = 0
    np.fill_diagonal(log_dist, 0)
    traj = md.load(topology)
    if np.any(input_contact_dist) == None and cutoff != np.inf:
        contact_dist, contact_pairs = md.compute_contacts(traj,
                                                       contacts = 'all',
                                                       scheme = 'closest-heavy',
                                                       ignore_nonprotein = False)

        squareform_dist = md.geometry.squareform(contact_dist, contact_pairs)[0]
        
    else:
        squareform_dist = input_contact_dist
        
    
    if cutoff == None:
        resultant_dist = np.multiply(log_dist, squareform_dist)
    elif cutoff == np.inf:
        resultant_dist=log_dist
    else:
        resultant_dist = np.multiply(log_dist, squareform_dist <= cutoff)
        
    return resultant_dist, log_dist, squareform_dist, correlation_array


def round_dict_values(input_dict):
    # Round all values in the dictionary to 3 decimal places
    rounded_dict = {k: round(v, 3) if isinstance(v, float) else v for k, v in input_dict.items()}
    return rounded_dict

def write_tcl_coloring(input_array, output_file = None):
    # Assuming the array is size n, where n is the number of residues to be colored

    input_lines = [] 

    input_lines.append('''# Assuming you're working with the first molecule loaded in VMD
set mol_id 0

# Get the total number of residues\n''')
            
    beta_line = 'set beta_values {'
    for value in input_array:
        beta_line += str(value)
        beta_line += ' '
    beta_line += '}\n'

    input_lines.append(beta_line)
    input_lines.append('# Get the total number of residues\n')
    input_lines.append('set num_residues {}\n'.format(len(input_array)))

    input_lines.append('''# Iterate over each residue
for {set i 0} {$i < $num_residues} {incr i} {
    # Get the beta value for this residue
    set beta_value [lindex $beta_values $i]
    set res [expr $i + 1]
    # Set the beta value for each atom in the residue
    # Assuming you're working with protein residues numbered sequentially
    set sel [atomselect $mol_id "resid $res"]
    $sel set beta $beta_value
    $sel delete
}
''')

    if output_file == None:
        print('# Note: No outputfile was provided, printing script in terminal ')
        for line in input_lines:
            print(line, end = '')
        # print... Write print code
    else: 
        try:
            f = open(output_file,'w+')
            for line in input_lines:
                f.write(line)
            f.close()
        except:
            print("Unable to write to file")

    return 

class CorrelationNetwork:
    def __init__(self):
        self.graph = nx.Graph()  
        # self.correlation_arrayin_correlation_array

        
    def graph_from_correlation(self, correlation_array, topology, cutoff = None, input_contact_dist = None):
        # Assuming adjacency_matrix is a square matrix (2D array or list of lists)
        adjacency_matrix, log_dist, squareform_dist, np_correlation_array = correlation_to_distance(correlation_array, topology, cutoff, input_contact_dist)
        self.graph = nx.from_numpy_array(adjacency_matrix)  # Add nodes
        self.adjacency= adjacency_matrix
        self.correlation = np_correlation_array
        self.log_dist = log_dist
        self.contact_squareform_dist = squareform_dist

        self.structure = md.load(topology)
        

                    
    def calculate_betweenness_centrality(self, weight='weight'):
        # Note: Ensure the 'weight' parameter matches how you've defined edge weights in the adjacency matrix
        self.betweenness_results = nx.betweenness_centrality(self.graph, weight=weight)
        return self.betweenness_results

    def manual_centrality(self):
        # Note: Likely makes sense not to use this manual function, although faster than the 
        # function in nx. That has been verified and will provide the correct units at least. 
        nodes_visted = []

        for i in tqdm(range(996)):
            d = nx.shortest_path(self.graph, target=i, weight='weight')
            # Extend the list with all values from the dictionary
            nodes_visted.extend(np.hstack(list(d.values())))

        # Concatenate all values into a NumPy array after the loop
        nodes_visted = np.array(nodes_visted)
        visted, node_counts = np.unique(nodes_visted, return_counts=True)
        self.node_counts = node_counts
        return node_counts  

    def shortest_path_to_residues(self, starting_residue, list_of_target_residues):
        # Subtract 1 to account for indexing of 0 based vs 1 based\
        targets = np.array(list_of_target_residues)-1
        source = starting_residue - 1 # ATP or dATP 
        print("Source   : {}".format(self.structure.top.residue(source)))

        pathways = []
        measured_lengths = []
        
        number_steps = []

        for resid in targets:
            
            print('From {} to target Resid: {}'.format(self.structure.top.residue(source), self.structure.top.residue(resid)))
            
            path = nx.shortest_path(self.graph, source = source, target = resid, weight = 'weight')
            number_steps.append(len(path))
            pathways.append(path)
            #print(ATP_path)
            for res in path:
                print(self.structure.top.residue(res), end = ' ')
            print()
            
            measured_lengths.append(nx.shortest_path_length(self.graph, source = source, target = resid, weight = 'weight'))
            print(measured_lengths[-1])

        return measured_lengths, number_steps, pathways
        

    

# def main():
        
#     example_top = '/home/marcus/Documents/SERCA/serca-MI/data/3W5A_ATP_rep2_local_5/3W5A_original_ATP_average.pdb'
#     example_correlation_array = '/home/marcus/Documents/SERCA/serca-MI/data/3W5A_ATP_rep2_local_5/3W5A_original_ATP_protein_align_0_6244.npy'
#     input_contact_dists = np.load('/home/marcus/Documents/SERCA/serca-MI/example_squareform_dist.npy')
#     #input_contact_dists = None 
#     example_network = CorrelationNetwork()
#     example_network.graph_from_correlation(example_correlation_array, example_top, input_contact_dist=input_contact_dists)
#     example_contact_dists = example_network.squareform_dist
#     # np.save('example_squareform_dist.npy',example_contact_dists)
#     #ATP_visits = example_network.manual_centrality()
#     write_tcl_coloring(example_network.correlation[995])


# if __name__ == "__main__":
#     main()
