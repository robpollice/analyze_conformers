# import packages
import pathlib as pl
import argparse as ap
import numpy as np
import matplotlib.pyplot as plt
import rdkit.Chem as rdc
import rdkit.Chem.rdMolTransforms as rdcmt

# define utility class
class ensemble:
    def __init__(self, path, name):
        # Initialize class variables
        self.name = name
        self.path = pl.Path(path)
        return
    
    def count_conformers(self):
        # Obtain number of conformers based on number of subdirectories
        self.conformers = sum(1 for pi in self.path.glob('*') if pi.is_dir())
        
        # Initiate array of principal moments for each conformer
        self.data = np.zeros((self.conformers,3))
        return
    
    def compute_principal_moments(self):
        # Compute principal moments for each conformer and store the results
        for ni in range(1, self.conformers + 1):
            # Open XYZ file
            fi = self.path / ni / 'conf.xyz'
            # Create mol object from XYZ
            mol = rdc.MolFromXYZFile(fi)
            # Retrieve conformer from mol object
            conformer = mol.GetConformer()
            # Compute principal moments
            _, principal_moments = rdcmt.ComputePrincipalAxesAndMoments(conformer, ignoreHs=False)
            self.data[ni-1, :] = principal_moments # Unit unclear
        return
    
    def generate_plots(self):
        # Plot scatter for 1st and 2nd principal moments of all conformers
        plt.scatter(self.data[:,0], self.data[:,1], c=list(range(self.conformers)), cmap='viridis_r')
        plt.xlabel('$1^{st}$ Principal Moment')
        plt.ylabel('$2^{nd}$ Principal Moment')
        cb = plt.colorbar()
        cb.ax.set_ylabel('Conformer Number')
        plt.grid()
        plot_name = 'principal_moments_1-2.pdf'
        plt.savefig(self.path / plot_name, format='pdf', bbox_inches='tight')
        plt.close()
        
        # Plot scatter for 1st and 3rd principal moments of all conformers
        plt.scatter(self.data[:,0], self.data[:,2], c=list(range(self.conformers)), cmap='viridis_r')
        plt.xlabel('$1^{st}$ Principal Moment')
        plt.ylabel('$3^{rd}$ Principal Moment')
        cb = plt.colorbar()
        cb.ax.set_ylabel('Conformer Number')
        plt.grid()
        plot_name = 'principal_moments_1-3.pdf'
        plt.savefig(self.path / plot_name, format='pdf', bbox_inches='tight')
        plt.close()
        
        # Plot scatter for 2nd and 3rd principal moments of all conformers
        plt.scatter(self.data[:,1], self.data[:,2], c=list(range(self.conformers)), cmap='viridis_r')
        plt.xlabel('$2^{nd}$ Principal Moment')
        plt.ylabel('$3^{rd}$ Principal Moment')
        cb = plt.colorbar()
        cb.ax.set_ylabel('Conformer Number')
        plt.grid()
        plot_name = 'principal_moments_2-3.pdf'
        plt.savefig(self.path / plot_name, format='pdf', bbox_inches='tight')
        plt.close()
        return

    def get_candidates(self, folded=True, percentile_delta=0.05):
        # Define parameters for conformer selection
        candidates = list()
        current_bound_value = 1 - float(folded)
        
        # Iterative candidate identification
        while len(candidates) == 0:
            current_bound_value += percentile_delta if folded == True else -1 * percentile_delta
            first_bound = np.quantile(self.data[:,0], 1 - current_bound_value)
            second_bound = np.quantile(self.data[:,1], current_bound_value)
            third_bound = np.quantile(self.data[:,2], current_bound_value)
            candidates = np.nonzero((self.data[:,0] < first_bound) & (self.data[:,1] > second_bound) & (self.data[:,2] > third_bound))[0]
        
        # Print results
        results_name = 'Folded' if folded == True else 'Unfolded'
        print(results_name + ' structures: ' + ', '.join(map(str, candidates + 1)))
        print(f'Best candidate: ' + str(candidates[0] + 1))
        print(f'Percentile used: {round(bound_value, 2)}')
        return
 

if __name__ == "__main__":
    # parse command line input
    parser = ap.ArgumentParser()
    parser.add_argument("path", type=str, help="Path containing crest results with individual conformers split into single files, each in a separate folder that is named with consecutive numbers.")
    parser.add_argument("--name", type=str, default="mol", help="Name of the molecule.")
    parser.add_argument("--plots", type=bool, default=True, help="Toggle whether scatter plots for principal moments of conformers should be generated.")
    args = parser.parse_args()
    
    # Create ensemble object
    conformer_ensemble = ensemble(args.path, args.name)
    
    # Process conformer ensemble
    conformer_ensemble.count_conformers()
    conformer_ensemble.compute_principal_moments()
    
    # Generate scatter plots of principal moments
    if args.plots == True:
        conformer_ensemble.generate_plots()
    
    # Get candidates for folded structures
    conformer_ensemble.get_candidates(folded=True, percentile_delta=0.05)
    
    # Get candidates for unfolded structures
    conformer_ensemble.get_candidates(folded=False, percentile_delta=0.10)
