import pandas as pd
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

# %% User inputs

# path to the file
file: str = r"C:\Users\AC03LH\OneDrive - Aalborg Universitet\Databases\VSA\Measurements\VON - HC Collet-08Dec2025-1,9997g.xls"
# file: str = 'VON - HFI1-test1-0,209g.xls'

# Dry weight of the sample in grams
dry_weight: float = 1.9997
# dry_weight: float = 0.209

# The VSA weighings differs a bit from the actual weight. Provide here a offset between the weights.
# An short analysis show that a general offset was around -0.0044 grams
vsa_offset: float = -0.0064  # TODO fix better corrigation


# %% Data treatement

def read(file: str) -> pd.DataFrame:
    df = pd.read_excel(file)
    return df


def init_data_treat(df: pd.DataFrame) -> pd.DataFrame:
    df['Weight (mg)'] = df['Weight (mg)'] - vsa_offset * 1000  # convert the offset in grams to an offset in miligrams
    df['% Moisture\nContent'] = (df[
                                     'Weight (mg)'] / 1000 - dry_weight) / dry_weight * 100  # Calculate the moisture content (%) in the sample
    return df


def DVS_handler(df: pd.DataFrame) -> list[pd.DataFrame]:
    df = df[df['Data\nType'] == 'Equilibration Point']
    return df


def DDI_handler(df: pd.DataFrame) -> list[pd.DataFrame]:
    df = df[df['Data\nType'] != 'Pre-Test']
    return df


def VSA_data_analyser(file: str, dry_weight: float, vsa_offset: float = 0) -> dict:
    data = read(file)
    data = init_data_treat(data)
    cycles = {}  # [stage][isotherm][data]
    print('Data contain these type of isotherms:')
    for stage in data['Stage'].unique():
        df_stage = data[data['Stage'] == stage].copy()
        method = df_stage.iloc[0]['Isotherm\nMethod']
        cycles[stage] = {}
        for isotherm, df_temp in df_stage.groupby('Sorption\nDirection'):
            print(f'Stage: {stage:<3} Method: {method:<5} Isotherm: {isotherm}')
            if method == 'DDI':
                cycles[stage][isotherm] = DDI_handler(df_temp)
            elif method == 'DVS':
                cycles[stage][isotherm] = DVS_handler(df_temp)
    return cycles


def VSA_plot(cycles: list, **kwargs) -> None:
    cycles = deepcopy(cycles)
    design = {1: {'Adsorption': {'c': 'blue', 'marker': '*', 'mec': 'black'},
                  'Desorption': {'c': 'blue', 'marker': '*', 'mec': 'black'}},
              2: {'Adsorption': {'c': 'orange', 'marker': '*', 'mec': 'black'},
                  'Desorption': {'c': 'orange', 'marker': '*', 'mec': 'black'}},
              3: {'Adsorption': {'c': 'tab:green', 'marker': '*', 'mec': 'black'},
                  'Desorption': {'c': 'tab:green', 'marker': '*', 'mec': 'black'}},
              4: {'Adsorption': {'c': 'red', 'marker': '*', 'mec': 'black'},
                  'Desorption': {'c': 'red', 'marker': '*', 'mec': 'black'}},
              5: {'Adsorption': {'c': 'brown', 'marker': '*', 'mec': 'black'},
                  'Desorption': {'c': 'brown', 'marker': '*', 'mec': 'black'}},
              6: {'Adsorption': {'c': 'cyan', 'marker': '*', 'mec': 'black'},
                  'Desorption': {'c': 'cyan', 'marker': '*', 'mec': 'black'}},
              }
    if not len(kwargs) == 0:
        msg = 'Selected plot options'
        print(f'{msg:-^30}')
    else:
        options = ['exclude']
        for key in kwargs.keys():
            if not key in options:
                raise Warning(f'{key} a valid option.')
    if 'exclude' in kwargs.keys():  # Exclude user specified stages
        n = kwargs['exclude']
        if type(n) == int:
            cycles.pop(n)
            print(f'Exclude stage {n:<2}')
        elif type(n) == list or type(n) == tuple:
            for i in n:
                if type(i) == int:
                    cycles.pop(i)
                    print(f'Exclude stage {i:<2}')
                elif type(i) == list or type(i) == tuple:
                    if not (len(i) == 2 and type(i[0]) == int and (i[1] == 'Adsorption' or i[1] == 'Desorption')):
                        raise Warning('The input for the exclude option is incorrect')
                    cycles[i[0]].pop(i[1])
                    print(f'Exclude stage {i[0]:<2} isotherm {i[1]:<11}')
    if 'layout' in kwargs.keys():  # Design options for the different stages
        for i in kwargs['layout']:
            print(f'Layout for stage {i} changed')
            design[i] = kwargs['layout'][i]

    fig = plt.figure(dpi=500, figsize=(7, 4))
    gs = fig.add_gridspec(ncols=1, nrows=1,
                          left=0.1, right=0.84, top=0.9, bottom=0.13)
    ax = fig.add_subplot(gs[0])
    for stage in cycles:
        for i, isotherm in enumerate(cycles[stage]):
            data_temp = cycles[stage][isotherm]
            x = data_temp['Water\nActivity']
            y = data_temp['% Moisture\nContent']
            if i == 0:
                ax.plot(x, y, label=f'Stage{stage}', **design[stage][isotherm])
            else:
                ax.plot(x, y, **design[stage][isotherm])
    ax.set_xticks(np.arange(0, 1.1, 0.1))
    ax.xaxis.set_major_formatter(PercentFormatter(1))
    ax.grid(color='black', alpha=0.3)
    fig.supxlabel('Water activity, [-]')
    fig.supylabel('% Moisture content by mass', x=0.01)
    # fig.suptitle(file)
    fig.legend(loc='center left', bbox_to_anchor=(0.84, 0.5))
    plt.show()


# %% Execute the scritp
if __name__ == '__main__':
    cycles = VSA_data_analyser(file, dry_weight, vsa_offset)

    exclude = (1, (3, 'Adsorption'))
    layout = {2: {'Adsorption': {'c': 'red', 'marker': '*', 'mec': 'black'},
                  'Desorption': {'c': 'tab:blue', 'marker': '*', 'mec': 'black'}}}
    exclude = (1,)
    VSA_plot(cycles, exclude=exclude, layout=layout)

    ad = []
    de = []
    for c in cycles:
        if c == 1:
            continue
        for i in cycles[c]:
            if i == 'Adsorption':
                ad.append(cycles[c][i][['Water\nActivity', '% Moisture\nContent']].groupby('Water\nActivity').first())
            else:
                de.append(cycles[c][i][['Water\nActivity', '% Moisture\nContent']].groupby('Water\nActivity').first())

ad = pd.concat(ad, axis=1).sort_index()
de = pd.concat(de, axis=1).sort_index()

# Bins og labels
bins = np.arange(0, 1.01, 0.01)
labels = [f'{0.01 * i:.2f}-{0.01 * i + 0.01:.2f}' for i in range(len(bins) - 1)]

# Kategorisering
ad_labels = pd.Series(pd.cut(ad.mean(axis=1).index, bins=bins, labels=labels, include_lowest=True).__array__(),
                      index=ad.index)
ad = pd.concat([ad.mean(axis=1), ad_labels], axis=1)
ad = ad.groupby(1).mean()
ad.rename(columns={0: 'ab'}, inplace=True)

de_labels = pd.Series(pd.cut(de.mean(axis=1).index, bins=bins, labels=labels, include_lowest=True).__array__(),
                      index=de.index)
de = pd.concat([de.mean(axis=1), de_labels], axis=1)
de = de.groupby(1).mean()
de.rename(columns={0: 'de'}, inplace=True)

final = pd.concat([ad, de], axis=1)
