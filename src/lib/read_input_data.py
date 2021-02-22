# Copyright (C) 2020-2021 Elisa Zuccolo, Eucentre Foundation
#
# HaselREC is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HaselREC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with HaselREC. If not, see <http://www.gnu.org/licenses/>.

def read_input_data(fileini):
    """
    """
    import sys
    import numpy as np

    input = {}
    with open(fileini) as fp:
        line = fp.readline()
        while line:
            if line.strip().find('=') >= 0:
                key, value = line.strip().split('=', 1)
                input[key.strip()] = value.strip()
            line = fp.readline()

    # %% Extract input parameters

    # Hazard parameters
    intensity_measures = [x.strip() for x in
                          input['intensity_measures'].strip('{}').split(',')]
    site_code = [x.strip() for x in input['site_code'].strip('{}').split(',')]
    site_code = np.array(site_code, dtype=int)
    rlz_code = [x.strip() for x in input['rlz_code'].strip('{}').split(',')]
    rlz_code = np.array(rlz_code, dtype=int)
    if len(rlz_code) != len(site_code):
        sys.exit(
            'Error: rlz_code must be an array of the same length of site_code')
    path_results_classical = input['path_results_classical']
    path_results_disagg = input['path_results_disagg']
    num_disagg = int(input['num_disagg'])
    num_classical = int(input['num_classical'])
    probability_of_exceedance_num = [x.strip() for x in input[
        'probability_of_exceedance_num'].strip('{}').split(',')]
    probability_of_exceedance_num = np.array(probability_of_exceedance_num,
                                             dtype=int)
    probability_of_exceedance = [x.strip() for x in
                                 input['probability_of_exceedance'].strip(
                                     '{}').split(',')]
    if (len(probability_of_exceedance_num) != len(
            probability_of_exceedance_num)):
        sys.exit(
            'Error: probability_of_exceedance_num must be of the same size of '
            'probability_of_exceedance')
    investigation_time = float(input['investigation_time'])

    # Conditional spectrum parameters
    target_periods = [x.strip() for x in
                      input['target_periods'].strip('[]').split(',')]
    target_periods = np.array(target_periods, dtype=float)

    tstar = np.zeros(len(intensity_measures))
    im_type = []
    im_type_lbl = []
    avg_periods = []

    for i in np.arange(len(intensity_measures)):
        if intensity_measures[i] == 'AvgSA':
            im_type.append('AvgSA')
            im_type_lbl.append(r'AvgSa')
            avg_periods = [x.strip() for x in
                           input['avg_periods'].strip('[]').split(',')]
            avg_periods = np.array(avg_periods, dtype=float)
        elif intensity_measures[i][0:2] == 'SA':
            im_type.append('SA')
            im_type_lbl.append(r'Sa(T)')
            tstar[i] = intensity_measures[i].strip('(,),SA')
            tstar[i] = float(tstar[i])
        elif intensity_measures[i][0:3] == 'PGA':
            im_type.append('PGA')
            im_type_lbl.append(r'PGA')
            tstar[i] = 0.0
        else:
            sys.exit('Error: this intensity measure type ' + str(
                intensity_measures[i]) + ' is not supported yet')

    # baker_jayaram or akkar
    corr_type = input['corr_type']
    gmpe_input = input['GMPE']  # array of GMPE according with sites?

    rake = float(input['rake'])

    vs30_input = [x.strip() for x in input['Vs30'].strip('{}').split(',')]
    if len(vs30_input) != len(site_code):
        sys.exit('Error: Vs30 must be an array of the same length of site_code')

    vs30type = [x.strip() for x in input['vs30Type'].strip('{}').split(',')]
    if len(vs30type) != len(site_code):
        sys.exit(
            'Error: vs30Type must be an array of the same length of site_code')

    hypo_depth = None
    try:
        hypo_depth = float(input['hypo_depth'])
    except KeyError:
        print(
            'Warning: if used, the hypocentral depth will be defined inside'
            ' the code')

    dip = None
    try:
        dip = float(input['dip'])
    except KeyError:
        print('Warning: if used, the dip angle will be defined inside the code')

    azimuth = None
    fhw = None
    try:
        azimuth = float(input['azimuth'])
    except KeyError:
        try:
            fhw = int(input['hanging_wall_flag'])
            if fhw != 1 and fhw != -1:
                sys.exit('Error: The hanging_wall_flag must be =1 or =-1')
        except KeyError:
            sys.exit(
                'Error: The azimuth or the hanging_wall_flag must be defined')

    z2pt5 = None
    try:
        z2pt5 = [x.strip() for x in input['z2pt5'].strip('{}').split(',')]
        if len(z2pt5) != len(site_code):
            sys.exit('Error: z2pt5 must be an array of the same length of'
                     ' site_code')
    except KeyError:
        print('Warning: if used, z2pt5 will be defined inside the code')

    z1pt0 = None
    try:
        z1pt0 = [x.strip() for x in input['z1pt0'].strip('{}').split(',')]
        if len(z1pt0) != len(site_code):
            sys.exit(
                'Error: z1pt0 must be an array of the same length of site_code')
    except KeyError:
        print('Warning: if used, z1pt0 will be defined inside the code')

    upper_sd = None
    try:
        upper_sd = float(input['upper_sd'])
    except KeyError:
        print('Warning: the upper_sd value will be defined inside the code')

    lower_sd = None
    try:
        lower_sd = float(input['lower_sd'])
    except KeyError:
        print('Warning: the lower_sd value will be defined inside the code')

    # Database parameters for screening recordings
    database_path = input['database_path']
    allowed_database = [x.strip() for x in
                        input['allowed_database'].strip('{}').split(',')]

    allowed_recs_vs30 = None
    try:
        # upper and lower bound of allowable Vs30 values
        allowed_recs_vs30 = [x.strip() for x in
                             input['allowedRecs_Vs30'].strip('[]').split(',')]
        allowed_recs_vs30 = np.array(allowed_recs_vs30, dtype=float)
    except KeyError:
        print('Warning: the vs30 range will be defined consistently withe the '
              'vs30 of the site')

    allowed_ec8_code = None
    try:
        allowed_ec8_code = [x.strip() for x in
                            input['allowedEC8code'].strip('{}').split(',')]
    except KeyError:
        print('Warning: the EC8 soil class will be defined consistently with'
              ' the EC8 soil class of the site ')

    try:
        maxsf_input = float(input['maxsf'])
    except ValueError:
        # The maximum allowable scale factor. They must be specified as a
        # function of probability_of_exceedance
        maxsf_input = [x.strip() for x in input['maxsf'].strip('{}').split(
            ',')]
        maxsf_input = np.array(maxsf_input, dtype=float)
        if len(probability_of_exceedance_num) != len(maxsf_input):
            sys.exit(
                'Error: maxsf must be of the same size of '
                'probability_of_exceedance')

    try:
        radius_dist_input = float(input['radius_dist'])
    except ValueError:
        radius_dist_input = [x.strip() for x in
                             input['radius_dist'].strip('{}').split(',')]  # km
        radius_dist_input = np.array(radius_dist_input, dtype=float)
        if len(probability_of_exceedance_num) != len(radius_dist_input):
            sys.exit('Error: radius_dist must be of the same size of '
                     'probability_of_exceedance')

    try:
        radius_mag_input = float(input['radius_mag'])
    except ValueError:
        radius_mag_input = [x.strip() for x in
                            input['radius_mag'].strip('{}').split(',')]
        radius_mag_input = np.array(radius_mag_input, dtype=float)
        if len(probability_of_exceedance_num) != len(radius_mag_input):
            sys.exit(
                'Error: radius_mag must be of the same size of '
                'probability_of_exceedance')

    # upper and lower bound of allowable depth values
    allowed_depth = [x.strip() for x in
                     input['allowed_depth'].strip('[]').split(
                         ',')]
    allowed_depth = np.array(allowed_depth, dtype=float)

    # Selection parameters
    # number of records to select ==> number of records to select since the code
    # search the database spectra most similar to each simulated spectrum
    n_gm = int(input['nGM'])
    random_seed = int(input['random_seed'])
    # number of iterations of the initial spectral simulation step to perform
    n_trials = int(input['nTrials'])
    # [Weights for error in mean, standard deviation and skewness] Used to find
    # the simulated spectra that best match the target from the statistically
    # simulated response spectra
    weights = [x.strip() for x in input['weights'].strip('{}').split(',')]
    weights = np.array(weights, dtype=float)
    # optimization parameters. Execution of incremental changes to the initially
    # selected ground motion set to further optimise its fit to the target
    # spectrum distribution.
    n_loop = int(input['nLoop'])  # Number of loops of optimization to perform
    # >0 to penalize selected spectra moire than 3 sigma from the target at any
    # period, =0 otherwise.
    penalty = float(input['penalty'])

    # Accelerogram folders
    path_nga_folder = input[
        'path_NGA_folder']  # NGA recordings have to be stored
    path_esm_folder = input['path_ESM_folder']
    # If not, found in the folder, ESM recording are authomatically downloaded
    # from internet, need to generate the file token.txt
    # At first you need to register at: https://tex.mi.ingv.it/
    # curl -X POST -F
    # 'message={"user_email": "email","user_password": "password"}'
    # "https://tex.mi.ingv.it/esmws/generate-signed-message/1/query" > token.txt

    # Output folder
    output_folder = input['output_folder']

    return (intensity_measures, site_code, rlz_code, path_results_classical,
            path_results_disagg, num_disagg, num_classical,
            probability_of_exceedance_num, probability_of_exceedance,
            investigation_time, target_periods, tstar, im_type, im_type_lbl,
            avg_periods, corr_type, gmpe_input, rake, vs30_input, vs30type,
            hypo_depth, dip, azimuth, fhw, z2pt5, z1pt0, upper_sd, lower_sd,
            database_path, allowed_database, allowed_recs_vs30,
            allowed_ec8_code, maxsf_input, radius_dist_input, 
            radius_mag_input, allowed_depth, n_gm, random_seed, n_trials,
            weights, n_loop, penalty, path_nga_folder, path_esm_folder,
            output_folder)
