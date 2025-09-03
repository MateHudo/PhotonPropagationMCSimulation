import numpy as np
import matplotlib.pyplot as plt

class LACLoader:
    def __init__(self, filename, material_density, material_name="MaterialName"):
        """
        Load mass attenuation coefficients from a txt file and convert to linear attenuation coefficients (LAC).
        Stores LACs for pair production (pp), Compton scattering (cs), photoelectric effect (pe), and their sum (total_simulation).

        Args:
            filename (str): Path to the MAC data file.
            material_density (float): Density of the material in g/cm^3.
            material_name (str): Name of the material for plot labeling.
        """
        data = np.loadtxt(filename)
        self.energy = data[:, 0]
        self.lac_cs = data[:, 2] * material_density
        self.lac_pe = data[:, 3] * material_density
        self.lac_pp = (data[:, 4] + data[:, 5]) * material_density
        self.lac_total_simulation = self.lac_cs + self.lac_pe + self.lac_pp
        self.material_name = material_name

    def get_lac_data(self):
        return self.energy, self.lac_pp, self.lac_cs, self.lac_pe, self.lac_total_simulation

    def plot(self, Emin=None, Emax=None):
        """
        Plot the LACs for the specified energy range.

        Args:
            Emin (float): Minimum energy to plot (default: min energy in data).
            Emax (float): Maximum energy to plot (default: max energy in data).
        """
        Emin = self.energy[0] if Emin is None else Emin
        Emax = self.energy[-1] if Emax is None else Emax
        mask = (self.energy >= Emin) & (self.energy <= Emax)
        plt.figure(figsize=(10, 6))
        plt.plot(self.energy[mask], self.lac_cs[mask], label='Compton Scattering (CS)', color='blue')
        plt.plot(self.energy[mask], self.lac_pe[mask], label='Photoelectric Effect (PE)', color='red')
        plt.plot(self.energy[mask], self.lac_pp[mask], label='Pair Production (PP)', color='green')
        plt.plot(self.energy[mask], self.lac_total_simulation[mask], label='Total LAC (CS + PE + PP)', color='black', linestyle='--')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('energija (MeV)')
        plt.ylabel('LAC (/cm)')
        plt.title(f'Linearni atenuacijski koeficienti za {self.material_name}')
        plt.legend()
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.tight_layout()
        plt.show()

    def plot_mean_free_path(self,Emin=None,Emax=None):
        """
        Plot the mean free path for the specified energy range.

        Args:
            Emin (float): Minimum energy to plot (default: min energy in data).
            Emax (float): Maximum energy to plot (default: max energy in data).
        """
        Emin = self.energy[0] if Emin is None else Emin
        Emax = self.energy[-1] if Emax is None else Emax
        mask = (self.energy >= Emin) & (self.energy <= Emax)
        mean_free_path = 1 / self.lac_total_simulation[mask]
        plt.figure(figsize=(10, 6))
        plt.plot(self.energy[mask], mean_free_path, label='Mean Free Path', color='purple')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('energija (MeV)')
        plt.ylabel('Srednja prosta pot (cm)')
        plt.title(f'Srednja prosta pot za {self.material_name}')
        plt.legend()
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.tight_layout()
        plt.show()

    def get_lac_at_energy(self, energy_value):
        """
        Get the LAC values at a specific energy using linear interpolation.

        Args:
            energy_value (float): Energy value in MeV.

        Returns:
            tuple: Interpolated LAC values (lac_pp, lac_cs, lac_pe, lac_total_simulation).
        """

        return np.interp(energy_value, self.energy, self.lac_total_simulation)
        