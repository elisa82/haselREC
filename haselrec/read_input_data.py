# Copyright (C) 2020-2021 Elisa Zuccolo, Eucentre Foundation
#
# haselREC is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# haselREC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with haselREC. If not, see <http://www.gnu.org/licenses/>.

def read_input_data(fileini):
    """
    Reads the input file for haselREC.
    An example is::

        [general]
        description = Selection for Italy - AvgSA
        selection_type = conditional-spectrum

        [hazard parameters]
        hazard_mode=0
        intensity_measures={AvgSA}
        site_code={0,1,2,3,4,5}
        rlz_code={0,0,0,0,0,0}
        probability_of_exceedance_num={2}
        probability_of_exceedance={0.1}
        path_results_classical=output_classical_AvgSA
        path_results_disagg=output_disagg_AvgSA
        num_disagg=107
        num_classical=108
        investigation_time=50

        [conditional spectrum parameters]
        target_periods=[0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9,1,2,3,4]
        corr_type=akkar
        GMPE=AkkarBommer2010
        avg_periods=[0.2,0.3,0.4,0.5,0.75,1.0]
        rake=0.
        Vs30={300,300,300,600,600,600}
        vs30Type={inferred,inferred,inferred,inferred,inferred,inferred}
        azimuth=50
        hypo_depth=10

        [database parameters for screening recordings]
        component=single-component
        database_path=../haselREC/GM-Records/
        allowed_database={ESM,NGA-West2,KiK-net}
        allowed_depth=[0,30]
        radius_dist={50}
        dist_range=[15,200]
        radius_mag={0.50}
        maxsf={3.0}

        [selection parameters]
        nGM=30
        nTrials = 10
        weights={0.6,0.3,0.1}
        nLoop=10
        penalty=10
        random_seed=333

        [accelerogram folders]
        path_NGA_folder=/rec/NGA2records
        path_ESM_folder=/rec/ESM
        path_kiknet_folder=/rec/KikNET

        [output folder]
        output_folder=output_acc_AvgSA

    The structure of the input file is similar to the structure of the input
    file for OpenQuake.

    **General - section**
    It contains the description of the run. It is not used by haselREC.
        - :code:`selection_type`: It can be conditional-spectrum or code-spectrum;

    **Hazard Parameters - section**
    It contains the hazard parameters, which must be defined according to the
    OpenQuake standards.

        - :code:`hazard_mode`: 0= hazard from OQ, 1= mean hazard given;
        - :code:`intensity_measures`: list of intensity measures to consider.
          They can be [`PGA`, `SA(T)` or `AvgSA`], where `T` is the spectral
          ordinate;
        - :code:`site_code`: list with site codes according to OpenQuake file
          names;
        - :code:`rlz_code`: list with realization codes according to OpenQuake
          file names (one for each site, :code:`site_code`);
        - :code:`probability_of_exceedance_num`: list with probability of
          exceedance numbers according to OpenQuake file names;
        - :code:`probability_of_exceedance`: list with the probability of
          exceedance associated to :code:`probability_of_exceedance_num`
          (only in case of :code:`hazard_mode`=0);
        - :code:`path_results_classical`: path to the folder containing the
          results of a classical PSHA performed by OpenQuake
          (only in case of :code:`hazard_mode`=0);
        - :code:`path_results_disagg`: path to the folder containing the
          results of a disaggregation analysis performed by OpenQuake
          (only in case of :code:`hazard_mode`=0);
        - :code:`num_disagg`: OpenQuake calculation ID containing
          the disaggregation results (only in case of :code:`hazard_mode`=0);
        - :code:`num_classical`: OpenQuake calculation ID containing
          the PSHA results (only in case of :code:`hazard_mode`=0);
        - :code:`investigation_time`: period of time (in years) used for PSHA in
          OpenQuake (only in case of :code:`hazard_mode`=0);
        - :code:`meanMag_disagg`: mean magnitude from disaggregation analysis 
          (only in case of :code:`hazard_mode`=1 or 
          :code:`selection_type`=code-spectrum. One value per site has to be specified);
        - :code:`meanDist_disagg`: mean distance from disaggregation analysis 
          (only in case of :code:`hazard_mode`=1 or;
          :code:`selection_type`=code-spectrum. One value per site has to be specified);
        - :code:`hazard_value`: hazard value from disaggregation analysis 
          (only in case of :code:`hazard_mode`=1)

    **Conditional Spectrum Parameters - section**

        - :code:`target_periods`: array of periods at which to compute the CS;
        - :code:`corr_type`: correlation relationship to be used for the
          computation of the CS. It can be [baker_jayaram or akkar];
        - :code:`GMPE`: name of the GMPE to be used for the the construction of
          the CS. It must be defined according to OpenQuake
          (see https://docs.openquake.org/oq-engine/master/openquake.hazardlib.
          gsim.html);
        - :code:`avg_periods`: range of periods to be used to compute AvgSA.
          It must be defined only when :code:`intensity_measures={AvgSA}`;
        - :code:`rake`: fault rake;
        - :code:`Vs30`: list of vs30 values (one for each target site,
          :code:`site_code`);
        - :code:`vs30Type`: list of vs30 sites (one for each target site,
          :code:`site_code`). It can be ["inferred" or "measured"];
        - :code:`azimuth`: (optional) source-to-site azimuth. It can be defined
          as an alternative to :code:`fhw`. If not defined, if required,
          it will be defined inside the code;
        - :code:`fhw`: (optional) hanging-wall flag. It can be [`=1` or `=-1`].
          If not defined, if required, it will be defined inside the code;
        - :code:`hypo_depth`: (optional) hypocentral depth. If not defined, if
          required, it will be defined inside the code;
        - :code:`dip`: (optional) fault dip. If not defined, if
          required, it will be defined inside the code;
        - :code:`z1pt0`: (optional) Depth to Vs=1.0 km/s. If not defined, if
          required, it will be defined inside the code;
        - :code:`z2pt5`: (optional) Depth to Vs=2.5 km/s. If not defined, if
          required, it will be defined inside the code;
        - :code:`upper_sd`: (optional) Upper seismogenic depth. If not defined,
          if required, it will be defined inside the code;
        - :code:`lower_sd`: (optional) Lower seismogenic depth. If not defined,
          if required, it will be defined inside the code;
        - :code:`allowed_recs_vs30`: (optional) Range of allowed vs30 values.
          If not defined, the vs30 range will be defined consistently withe the
          vs30 of the site;
        - :code:`allowed_ec8_code`: (optional). List of allowed EC8 soil class
          codes. It can be ['A', 'B', 'C', 'D', 'E' or 'All']. If not defined,
          the EC8 soil class will be defined consistently withe the EC8 soil
          class of the site;

    **Code Spectrum Parameters - section**
        - :code:`code_spectrum_file`: path to the files containing the code spectrum (in g).
          The number of files should be equal to the number of sites
        - :code:`period_range_spectrumcompatibility`: period interval for 
          spectrum compatibility (sec)
        - :code:`threshold_up`: maximum threshold for the maximum positive 
          difference (%)
        - :code:`threshold_low`: = minimum tthreshold for the maximum negative
          dfference (%)
        - :code:`allowed_recs_vs30`: (optional) Range of allowed vs30 values.
          If not defined, the vs30 range will be defined consistently withe the
          vs30 of the site;
        - :code:`allowed_ec8_code`: (optional). List of allowed EC8 soil class
          codes. It can be ['A', 'B', 'C', 'D', 'E' or 'All']. If not defined,
          the EC8 soil class will be defined consistently withe the EC8 soil
          class of the site;

    **Database Parameters For Screening Recordings - section**

        - :code:`component`: it can be single-component or two-component
        - :code:`database_path`: path to the folder containing the strong motion
          database;
        - :code:`allowed_database`: list of databases to consider for record
          selection. They can be ['NGA-West2' or 'ESM' or 'KiK-net'];
        - :code:`allowed_depth`: upper and lower bound of allowable depths;
        - :code:`radius_dist`: list of radius values to be used for the
          definition of the allowable distance values. They must be
          specified for each probability of exceedance (:code:
          `probability_of_exceedance`);
        - :code:`dist_range`: minimum and maximum distance values to be used
          when defining the distance range according to (:code:`radius_dist`)
        - :code:`radius_mag`: list of radius values to be used for the
          definition of the allowable magnitude values. They must be
          specified for each probability of exceedance
          (:code:`probability_of_exceedance`);
        - :code:`radius_mag_type`: list of radius types to be used for the
          definition of the allowable magnitude values. They must be
          specified for each probability of exceedance. They can be "both",
          "left", "right" and they indicate if the radius should be applied
          in both directions or only in one direction.
        - :code:`radius_dist_type`: list of radius types to be used for the
          definition of the allowable magnitude values. They must be
          specified for each probability of exceedance. They can be "both",
          "left", "right" and they indicate if the radius should be applied
          in both directions or only in one direction.
        - :code:`maxsf`: list of maximum allowable scale factor. They must be
          specified for each probability of exceedance (:code:`probability_of_exceedance`);
        - :code:`correlated_motion`: it can be yes or no. If no, 
          exclude ground motions from the same stations and earthquake. The 
          default is yes (i.e. consider also recordings from same stations and 
          earthquakes)

    **Selection Parameters - section**

        - :code:`nGM`: number of records to select;
        - :code:`nTrials`: number of iterations of the initial spectral
          simulation step to perform;
        - :code:`weights`: {weight for error in mean, weight for error in
          standard deviation, weight for error in skewness};
        - :code:`nLoop`: number of loops of optimization to perform;
        - :code:`penalty`: >0 to penalize selected spectra more than 3 sigma
          from the target at any period, =0 otherwise;
        - :code:`random_seed`: random seed number to simulate response spectra
          for initial matching;
        - :code:`scaling_vertical_component`: decide if the vertical component has to 
          be scaled or not. It can be "yes" or "no"

    **Accelerogram Folders - section**

        - :code:`path_NGA_folder`: folder with NGA-Wes2 recordings;
        - :code:`path_ESM_folder`: folder with ESM recordings (it can be empty);
        - :code:`path_kiknet_folder`: folder with KiK-net recordings (it can be empty);

    **Output Folder - section**

        - :code:`output_folder`: path to the output folder.
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

    selection_type = input['selection_type']

    # %% Extract input parameters

    # Hazard parameters

    site_code = [x.strip() for x in input['site_code'].strip('{}').split(',')]
    site_code = np.array(site_code, dtype=int)
    probability_of_exceedance_num = [x.strip() for x in input[
        'probability_of_exceedance_num'].strip('{}').split(',')]
    probability_of_exceedance_num = np.array(probability_of_exceedance_num,
                                             dtype=int)

    if( selection_type == 'conditional-spectrum'): 
        hazard_mode=int(input['hazard_mode'])
        intensity_measures = [x.strip() for x in
                              input['intensity_measures'].strip('{}').split(',')]
        rlz_code = [x.strip() for x in input['rlz_code'].strip('{}').split(',')]
        rlz_code = np.array(rlz_code, dtype=int)
        if len(rlz_code) != len(site_code):
            sys.exit(
                'Error: rlz_code must be an array of the same length of site_code')
        if(hazard_mode==0):
            path_results_classical = input['path_results_classical']
            path_results_disagg = input['path_results_disagg']
            num_disagg = int(input['num_disagg'])
            num_classical = int(input['num_classical'])
            probability_of_exceedance = [x.strip() for x in
                                         input['probability_of_exceedance'].strip(
                                             '{}').split(',')]
            if (len(probability_of_exceedance_num) != len(
                    probability_of_exceedance_num)):
                sys.exit(
                    'Error: probability_of_exceedance_num must be of the same size of '
                    'probability_of_exceedance')
            investigation_time = float(input['investigation_time'])
            meanMag_disagg=[]
            meanDist_disagg=[]
            hazard_value=[]
        else:
            meanMag_disagg = [x.strip() for x in input['meanMag_disagg'].strip('{}').split(',')]
            meanMag_disagg = np.array(meanMag_disagg, dtype=float)
            meanDist_disagg = [x.strip() for x in input['meanDist_disagg'].strip('{}').split(',')]
            meanDist_disagg = np.array(meanDist_disagg, dtype=float)
            hazard_value=float(input['hazard_value'])
            path_results_classical = []
            path_results_disagg = []
            num_disagg = []
            num_classical = []
            probability_of_exceedance = []
            investigation_time=[]
    elif(selection_type == 'code-spectrum'): 
        meanMag_disagg = [x.strip() for x in input['meanMag_disagg'].strip('{}').split(',')]
        meanMag_disagg = np.array(meanMag_disagg, dtype=float)
        meanDist_disagg = [x.strip() for x in input['meanDist_disagg'].strip('{}').split(',')]
        meanDist_disagg = np.array(meanDist_disagg, dtype=float)
        hazard_mode=[]
        intensity_measures = [0.0]
        hazard_value=[]
        path_results_classical = []
        path_results_disagg = []
        num_disagg = []
        num_classical = []
        probability_of_exceedance = [0]
        probability_of_exceedance = np.array(probability_of_exceedance)
        investigation_time=[]
        rlz_code=np.zeros((len(site_code)))
    else:
        sys.exit('Error: this selection type ' + str(selection_type) 
                + ' is not supported')

    #Spectrum parameters 

    target_periods = []
    im_type = []
    im_type_lbl = []
    avg_periods = []
    rake = None
    gmpe_input = None
    corr_type = None
    vs30_input = []
    vs30type = []
    hypo_depth = None
    dip = None
    azimuth = None
    fhw = None
    z2pt5 = None
    z1pt0 = None
    upper_sd = None
    lower_sd = None
    threshold_up = None
    threshold_low = None
    code_spectrum_file = None
    period_range_spectrumcompatibility = []
    tstar = np.zeros(len(intensity_measures))

    if( selection_type == 'conditional-spectrum'):

    # Conditional spectrum parameters
        target_periods = [x.strip() for x in
                          input['target_periods'].strip('[]').split(',')]
        target_periods = np.array(target_periods, dtype=float)

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
        gmpe_input = input['GMPE']
        rake = float(input['rake'])

        vs30_input = [x.strip() for x in input['Vs30'].strip('{}').split(',')]
        if len(vs30_input) != len(site_code):
            sys.exit('Error: Vs30 must be an array of the same length of site_code')

        vs30type = [x.strip() for x in input['vs30Type'].strip('{}').split(',')]
        if len(vs30type) != len(site_code):
            sys.exit(
                'Error: vs30Type must be an array of the same length of site_code')

        try:
            hypo_depth = float(input['hypo_depth'])
        except KeyError:
            print(
                'Warning: if used, the hypocentral depth will be defined inside'
                ' the code')

        try:
            dip = float(input['dip'])
        except KeyError:
            print('Warning: if used, the dip angle will be defined inside the code')

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

        try:
            z2pt5 = [x.strip() for x in input['z2pt5'].strip('{}').split(',')]
            if len(z2pt5) != len(site_code):
                sys.exit('Error: z2pt5 must be an array of the same length of'
                         ' site_code')
        except KeyError:
            print('Warning: if used, z2pt5 will be defined inside the code')

        try:
            z1pt0 = [x.strip() for x in input['z1pt0'].strip('{}').split(',')]
            if len(z1pt0) != len(site_code):
                sys.exit(
                    'Error: z1pt0 must be an array of the same length of site_code')
        except KeyError:
            print('Warning: if used, z1pt0 will be defined inside the code')

        try:
            upper_sd = float(input['upper_sd'])
        except KeyError:
            print('Warning: the upper_sd value will be defined inside the code')

        try:
            lower_sd = float(input['lower_sd'])
        except KeyError:
            print('Warning: the lower_sd value will be defined inside the code')

        tstar1=-1
        tstar2=-1

    # Code spectrum parameters
    if( selection_type == 'code-spectrum'):
        scaling_period=input['scaling_period']
        try:
            tstar[0]=float(scaling_period)
            tstar1=-1
            tstar2=-1
        except ValueError:
            tstar_temp=[x.strip() for x in
                             scaling_period.strip('[]').split(',')]
            tstar1=float(tstar_temp[0])
            tstar2=float(tstar_temp[1])

        im_type_lbl.append(r'None')
        code_spectrum_file = [x.strip() for x in input['code_spectrum_file'].strip('{}').split(',')]
        if len(code_spectrum_file) != len(site_code):
            sys.exit(
                'Error: code_spectrum_file must be an array of the same length of site_code')
        period_range_spectrumcompatibility = [x.strip() for x in
            input['period_range_spectrumcompatibility'].strip('[]').split(
                         ',')]
        period_range_spectrumcompatibility = np.array(period_range_spectrumcompatibility, dtype=float)
        threshold_up = float(input['threshold_up'])
        threshold_low = float(input['threshold_low'])
        
    # Database parameters for screening recordings
    component = input['component']
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

    correlated_motion='yes'
    try:
        correlated_motion = input['correlated_motion']
    except KeyError:
         print('The selection could include recordings associated to not '
                 'unique stations and earthquakes')

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
        radius_dist_type_input = input['radius_dist_type']
    except ValueError:
        radius_dist_type_input = [x.strip() for x in
                             input['radius_dist_type'].strip('{}').split(',')]
        if len(probability_of_exceedance_num) != len(radius_dist_input):
                sys.exit('Error: radius_dist_type must be of the same size of '
                         'radius_dist')

    dist_range_input = None
    try:
        #Distance range allowed
        dist_range_input = [x.strip() for x in
                             input['dist_range'].strip('[]').split(',')]
        dist_range_input = np.array(dist_range_input, dtype=float)
    except KeyError:
        print('Warning: an ad-hoc dist range is not defined by the user. Therefore no additional contraints on the distance range will be considered. Only radius_dist will be adopted')

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

    try:
        radius_mag_type_input = input['radius_mag_type']
    except ValueError:
        radius_mag_type_input = [x.strip() for x in
                             input['radius_mag_type'].strip('{}').split(',')]
        if len(probability_of_exceedance_num) != len(radius_mag_input):
            sys.exit('Error: radius_mag_type must be of the same size of '
                    'radius_mag')

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

    try:
        vertical_component = str(input['scaling_vertical_component'])
    except ValueError:
        vertical_component = 'no' 

    # Accelerogram folders
    path_nga_folder = input['path_NGA_folder']# NGA recordings have to be stored
    path_esm_folder = input['path_ESM_folder']
    path_kiknet_folder = input['path_kiknet_folder']

    # Output folder
    output_folder = input['output_folder']

    return (intensity_measures, site_code, rlz_code, path_results_classical,
            path_results_disagg, num_disagg, num_classical,
            probability_of_exceedance_num, probability_of_exceedance,
            investigation_time, target_periods, tstar, im_type, im_type_lbl,
            avg_periods, corr_type, gmpe_input, rake, vs30_input, vs30type,
            hypo_depth, dip, azimuth, fhw, z2pt5, z1pt0, upper_sd, lower_sd,
            database_path, allowed_database, allowed_recs_vs30,
            allowed_ec8_code, maxsf_input, radius_dist_input, dist_range_input,
            radius_mag_input, allowed_depth, n_gm, random_seed, n_trials,
            weights, n_loop, penalty, path_nga_folder, path_esm_folder,
            output_folder, meanMag_disagg, meanDist_disagg, hazard_value, 
            hazard_mode, component, correlated_motion, code_spectrum_file,
            period_range_spectrumcompatibility, threshold_up, threshold_low,
            selection_type, path_kiknet_folder, 
            radius_dist_type_input, radius_mag_type_input, vertical_component,
            tstar1,tstar2)

