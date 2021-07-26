"""
Generator of synthetic microbial culture data under anherobic, chemotrophic
growth. Growth and yield are controlled by user-provided hyperparameters.

Antonio Paiva, CSR
"""

import numpy as np
import matplotlib.pyplot as plt
import pdb
from pprint import pprint


class SRBGrowthSim():
    """SRB Microbial culture chemotrophic growth data generator.

    Assumed "experimental" conditions:
        Chemotrophic growth of SRB in 3mL solution, with lactate acting
        both the as electron donor and carbon source.

        Catabolic reaction:
            -[lactate] - 0.5[sulfate] + [acetate] + [bicarbonate]
                + 0.5[dihydrogen sulfide] = 0
        Anabolic reaction:
            -0.35[lactate] - 0.2[ammonia] - 0.1[proton] + [biomass]
                + 0.05[acetate] + 0.4[water] = 0

        The biomass yield is the "efficiency" with which the electron donor
        (lactate in this case) is converted into biomass. This determines
        the multiplier onto the catabolic reaction equation for combining
        these two reactions into the overall macrochemical reaction equation.

        By stoichiometric balance, the biomass yield is given by:
            Y_DX = 1/(gamma + 0.35)
        where gamma is the catabolic multipler and 0.35 the stoichiometric
        coefficient of lactate in the anabolic reaction. Hence,
            gamma = (1/Y_DX) - 0.35

    Arguments:
        biomass_yield_curve: function
            Function the returns the biomass yield depending on temperature
            and starting lactate concentration.

        cell_weight_curve: function
            Function that yields the cell weight depending on the temperature
            and starting lactate concentration.

        concentrations_noise: dict
            Standard deviation of normal distributed measurement noise added
            to the observed concentrations of the species given as keys of
            the dict. Since the concentrations would be measured at the
            beginnning and end of the experiments, this is added with a
            `sqrt(2)` scaling to the delta in concentrations.

        cell_count_noise: float
            Standard deviation of normal distributed cell count error.
            This is assumed to be somewhat proportional to the number of
            cells, and thus determined by the final cell count
            (i.e., no `sqrt(2)` scaling).
    """

    def __init__(self,
                 biomass_yield_curve,
                 cell_weight_curve,
                 concentrations_noise,
                 cell_count_noise):

        self.biomass_yield_curve = biomass_yield_curve
        self.cell_weight_curve = cell_weight_curve
        self.concentrations_noise = concentrations_noise
        self.cell_count_noise = cell_count_noise

        # general assumptions / simulation parameters
        self.biomass_molar_mass = 24.62     # g/C-mol
        self.lactate_final = 0.01    # amount of residual lactate
        self.volume = 3e-3           # volume of the vial/solution, in liters
        #self.initial_cell_weight = 1.8e-13  # in grams

    def generate(self, temperature, initial_cell_count,
                 initial_lactate_concentration=30.,
                 initial_sulfate_concentration=30.,
                 initial_acetate_concentration=0.1,
                 initial_bicarbonate_concentration=20.,
                 initial_sulfide_concentration=0.2):
        """Generate growth data for "one experiment".

        Returns:
            1. Latent parameters (initial and final biomass, lactate consumed
               in catabolic & anabolic reactions)
            2. Initial conditions (cell counts and concentrations of lactate,
               sulfate, acetate, bicarbonate, and sulfide)
            3. Change in conditions (negative [resp. positive] numbers mean
               that the compound was consumed [resp. produced]
        """

        # sample initial concentrations by adding noise
        # add noise to lactate now but not in the final reading
        initial_conditions = {
            'n_cells': initial_cell_count,
            'lactate': initial_lactate_concentration,
            'sulfate': initial_sulfate_concentration,
            'acetate': initial_acetate_concentration,
            'sulfide': initial_sulfide_concentration,
            'bicarbonate': initial_bicarbonate_concentration}

        # determine general growth conditions according to the
        # user-provided relationships
        yield_ = self.biomass_yield_curve(temperature,
                                          initial_conditions['lactate'])
        cell_weight = self.cell_weight_curve(temperature,
                                             initial_conditions['lactate'])
        gamma = self.gamma(yield_)
        #print('DEBUG: yield =', yield_)
        #print('DEBUG: gamma =', gamma)

        # determine how much biomass was produced for the target yield and
        # how much lactate was used in catabolic and anabolic reactions
        # **units: millimoles**
        lactate_consumed = (initial_lactate_concentration
                            - self.lactate_final) * (self.volume/1e3)
        total_biomass_produced = yield_ * lactate_consumed
        lactate_ana = 0.35 * total_biomass_produced
        lactate_cat = lactate_consumed - lactate_ana

        lactate_cat_conc = 1e3*(lactate_cat / self.volume)  # concentration [mM]
        lactate_ana_conc = 1e3*(lactate_ana / self.volume)  # concentration [mM]

        # determine how much of the biomass production is observed in the
        # cell counts vs changed in the cell weights
        initial_biomass = initial_cell_count * cell_weight \
                            / self.biomass_molar_mass  # final units: C-moles
        delta_cell_count = total_biomass_produced \
            * (self.biomass_molar_mass / cell_weight)

        # determine final conditions
        delta_conditions = {
            'n_cells': self.add_delta_noise(delta_cell_count, 'n_cells'),
            'lactate': self.add_delta_noise(-lactate_cat_conc-lactate_ana_conc, 'lactate'),
            'sulfate': self.add_delta_noise(-0.5*lactate_cat_conc, 'sulfate'),
            'acetate': self.add_delta_noise(lactate_cat_conc, 'acetate'),
            'sulfide': self.add_delta_noise(0.5*lactate_cat_conc, 'sulfide'),
            'bicarbonate': self.add_delta_noise(
                lactate_cat_conc + (0.05/0.35)*lactate_ana_conc, 'bicarbonate')}

        # first return arg are latent experiment parameters
        return ((initial_biomass, initial_biomass+total_biomass_produced,
                 -lactate_cat_conc, -lactate_ana_conc, gamma),
                initial_conditions, delta_conditions)

    def add_delta_noise(self, x, species):
        if species == 'n_cells':
            n  = x + np.round(self.cell_count_noise * np.random.randn())
            return n if n > 0 else 0
        elif species in self.concentrations_noise:
            return x + np.sqrt(2)*self.concentrations_noise[species] \
                        * np.random.randn()
        else:
            warning(f'Did not find noise stddev specification for: {species}')
            return x

    def gamma(self, yield_):
        return (1/yield_) - 0.35


if __name__ == '__main__':
    # example code on how to use the simulator

    # auxiliary functions for specifying the curves
    s = lambda x: 1./(1 + np.exp(-x))  # sigmoid
    g = lambda x, e=2: np.exp(-0.5*(np.abs(x)**e))  # gaussian

    # generation curves
    def growth_rate_function(t, l0):
        return np.maximum(1e-6, 0.021*g((t-18)/12, e=1.6)
                                     *s(-1.1*(t-18-8/1.1)) - 0.001)

    def biomass_yield_function(t, l0):
        return np.maximum(1e-6, 0.24 - 0.008*t)

    def cell_weight_function(t, l0):
        #p = np.polyfit([0, 20, 32, 35, 40], [2.5, 1.6, 1., 1.5, 2.5], deg=3)
        #return 1.8e-13*np.polyval(p, t)
        return 1.8e-13*(t > -10)

    if True:
        t = np.linspace(-2, 28, 81)
        fig, ax = plt.subplots(ncols=3, figsize=(9,3))
        ax[0].plot(t, growth_rate_function(t, None))
        ax[0].set_title('growth_rate_function')
        ax[1].plot(t, biomass_yield_function(t, None))
        ax[1].set_title('biomass_yield_function')
        ax[2].plot(t, cell_weight_function(t, None))
        ax[2].set_title('cell_weight_function')
        plt.tight_layout()
        plt.savefig('curves.png')
        plt.show(block=False)

    # initialize the data generator
    sim = SRBGrowthSim(
                 biomass_yield_curve=biomass_yield_function,
                 cell_weight_curve=cell_weight_function,
                 concentrations_noise={
                     'lactate': 0.3,
                     'sulfate': 0.6,
                     'acetate': 0.4,
                     'bicarbonate': 2.0,
                     'sulfide': 0.6},
                 cell_count_noise=1e7)

    latent, initial_conditions, delta_conditions = sim.generate(20, 1e7)

    print('latent =', latent)
    print('initial_conditions ='); pprint(initial_conditions)
    print('delta_conditions ='); pprint(delta_conditions)
