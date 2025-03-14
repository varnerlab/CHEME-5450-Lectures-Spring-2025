# -- PRIVATE HELPER METHODS ------------------------------------------------------------------------------- #
function compute_transcription_rate(parameter_dictionary::Dict{Symbol,Any})

    # compute rX_hat bound -
    VX = parameter_dictionary[:VX]
    length_factor_transcription = parameter_dictionary[:length_factor_transcription]
    gene_concentration = parameter_dictionary[:gene_concentration]
    X_tau_factor = parameter_dictionary[:X_tau_factor]
    KX = parameter_dictionary[:KX]


    # compute the control variable value -
    induction_parameter_dictionary = parameter_dictionary[:induction_parameter_dictionary]
    K1 = induction_parameter_dictionary[:K1]
    K2 = induction_parameter_dictionary[:K2]
    n = induction_parameter_dictionary[:coop_parameter]
    k_binding = induction_parameter_dictionary[:k_binding]
    I = induction_parameter_dictionary[:inducer_concentration]
    f_function = (I^n)/(k_binding^n+I^n)
    u_variable = (K1+K2*f_function)/(1+K1+K2*f_function)


    # kinetic limit -
    transcription_rate = VX*length_factor_transcription*(gene_concentration/(KX*X_tau_factor+(1+X_tau_factor)*gene_concentration));

    # return transcription rate -
    return transcription_rate*u_variable
end

function compute_translation_rate(rX_hat, parameter_dictionary::Dict{Symbol,Any})

    # Get some stuff from the parameter dictionary -
    VL = parameter_dictionary[:VL]
    length_factor_transcription = parameter_dictionary[:length_factor_translation]
    L_tau_factor = parameter_dictionary[:L_tau_factor]
    KL = parameter_dictionary[:KL]

    # compute the steady-state mRNA level -
    kdX = parameter_dictionary[:kdX]
    m_ss = (rX_hat/kdX)

    # compute the rL_hat -
    rL_hat = VL*length_factor_transcription*(m_ss/(KL*L_tau_factor+(1+L_tau_factor)*m_ss))

    # return -
    return rL_hat
end

function open_default_bounds_array(parameter_dictionary::Dict{Symbol,Any})

    # get the default bounds array -
    flux_bounds_array = deepcopy(parameter_dictionary[:flux_bounds_array])

    # ok, so the bounds that matters are the X (R2) and L (R5) bounds -
    rX_hat = compute_transcription_rate(parameter_dictionary)
    rL_hat = compute_translation_rate(rX_hat,parameter_dictionary)

    # update the bounds -
    # transcription rate (eqaulity constrained to rX_hat)
    flux_bounds_array[2,1] = 0          # lower bound on transcription
    flux_bounds_array[2,2] = rX_hat     # upper bound on transcription

    # translation rate (0 -> rL_hat)
    flux_bounds_array[5,1] = 0          # lower bound is no translation
    flux_bounds_array[5,2] = rL_hat     # upper bound is rL_hat

    # package -
    parameter_dictionary[:flux_bounds_array] = flux_bounds_array

    # return -
    return parameter_dictionary

end

function update_default_bounds_array(parameter_dictionary::Dict{Symbol,Any})

    # get the default bounds array -
    flux_bounds_array = deepcopy(parameter_dictionary[:flux_bounds_array])

    # ok, so the bounds that matters are the X (R2) and L (R5) bounds -
    rX_hat = compute_transcription_rate(parameter_dictionary)
    rL_hat = compute_translation_rate(rX_hat,parameter_dictionary)

    # update the bounds -
    # transcription rate (eqaulity constrained to rX_hat)
    flux_bounds_array[2,1] = rX_hat     # lower bound on transcription
    flux_bounds_array[2,2] = rX_hat     # upper bound on transcription

    # translation rate (0 -> rL_hat)
    flux_bounds_array[5,1] = 0          # lower bound is no translation
    flux_bounds_array[5,2] = rL_hat     # upper bound is rL_hat

    # package -
    parameter_dictionary[:flux_bounds_array] = flux_bounds_array

    # return -
    return parameter_dictionary
end

function generate_parameter_dictionary(path_to_parameter_file::String)

    # some constants -
    av_number = 6.02e23

    # load the parameter JSON file -
    simulation_json_tree = JSON.parsefile(path_to_parameter_file)

    # initialize -
    parameter_dictionary = Dict{Symbol,Any}()

    # do some stuff, and populate the dictionary -
    # stoichiometrix_matrix = readdlm(simulation_json_tree["paths_to_model_files"]["path_to_stochiometric_matrix_file"])
    # flux_bounds_array = readdlm(simulation_json_tree["paths_to_model_files"]["path_to_flux_bounds_file"])

    # # remove the unbalanced species ([e]'s start at 18, see Debug.txt')
    # balanced_stoichiometrix_matrix = stoichiometrix_matrix[1:17,:]

    # # get system size so I can make the species bounds array -
    # (number_of_species, number_of_fluxes) = size(balanced_stoichiometrix_matrix)
    # species_bounds_array = zeros(number_of_species,2)

    # # by default, objective_coefficient_array is zeros -
    # objective_coefficient_array = zeros(number_of_fluxes)

    # calculate the gene concentration -
    gene_copy_number = parse(Float64,simulation_json_tree["simulation_paramaters"]["gene_copy_number"])
    volume_in_L = (1/1e6)*parse(Float64,simulation_json_tree["simulation_paramaters"]["cf_reaction_volume_in_muL"])
    gene_concentration = (gene_copy_number)*(1/volume_in_L)*(1/av_number)*(1e6) # μM
    gene_concentration = 5.0*(1e6/1e9) # μM

    # calculate the length factors -
    length_of_gene_in_nt = parse(Float64,simulation_json_tree["simulation_paramaters"]["length_of_gene_in_nt"])
    characteristic_gene_length = parse(Float64,simulation_json_tree["biophysical_parameters"]["average_transcript_length"]["value"])
    length_factor_transcription = (characteristic_gene_length/length_of_gene_in_nt)

    length_of_prot_in_aa = parse(Float64,simulation_json_tree["simulation_paramaters"]["length_of_protein_in_aa"])
    characteristic_prot_length = parse(Float64,simulation_json_tree["biophysical_parameters"]["average_protein_length"]["value"])
    length_factor_translation = (characteristic_prot_length/length_of_prot_in_aa)

    # calculate the concentrations of RNAP (RX) and Ribosome (RT) -
    cf_dilution_factor = parse(Float64,simulation_json_tree["simulation_paramaters"]["cf_dilution_factor"])
    volume_of_single_cell = parse(Float64,simulation_json_tree["biophysical_parameters"]["volume_of_single_cell"]["value"]) # L

    RX_copy_number = parse(Float64,simulation_json_tree["biophysical_parameters"]["copies_of_rnapII_per_cell"]["value"]) # #/cell
    RX_concentration = (RX_copy_number)*(1/av_number)*(1/volume_of_single_cell)*(1/cf_dilution_factor)*(1e6) # μM

    RL_copy_number = parse(Float64,simulation_json_tree["biophysical_parameters"]["copies_of_ribosome_per_cell"]["value"]) # #/cell
    RL_concentration = (RL_copy_number)*(1/av_number)*(1/volume_of_single_cell)*(1/cf_dilution_factor)*(1e6) # μM

    # calculate the time time constant for transcription -
    max_transcription_rate = parse(Float64,simulation_json_tree["biophysical_parameters"]["transcription_elongation_rate"]["value"])
    transcription_initiation_time_contstant = parse(Float64,simulation_json_tree["biophysical_parameters"]["transcription_initiation_time_constant"]["value"])
    kcat_transcription = max_transcription_rate*(3600/length_of_gene_in_nt)                         # hr^-1
    kcat_transcription_initiation = (1/transcription_initiation_time_contstant)*(3600)              # hr^-1
    X_tau_factor = (kcat_transcription)/(kcat_transcription_initiation)                             # dimensionless

    # calculate the time constant for translation -
    max_translation_rate = parse(Float64,simulation_json_tree["biophysical_parameters"]["translation_elongation_rate"]["value"])
    translation_initiation_time_contstant = parse(Float64,simulation_json_tree["biophysical_parameters"]["translation_initiation_time_constant"]["value"])
    kcat_translation = max_translation_rate*(3600/length_of_prot_in_aa)                             # hr^-1
    kcat_translation_initiation = (1/translation_initiation_time_contstant)*(3600)                  # hr^-1
    L_tau_factor = (kcat_translation)/(kcat_translation_initiation)                                 # dimensionless

    # calculate the degradation constants for mRNA and protein -
    mRNA_half_life = parse(Float64,simulation_json_tree["biophysical_parameters"]["average_mRNA_half_life"]["value"]) # hr
    kdX = -(1/mRNA_half_life)*log(0.5)

    protein_half_life = parse(Float64,simulation_json_tree["biophysical_parameters"]["average_protein_half_life"]["value"]) # hr
    kdL = -(1/protein_half_life)*log(0.5)

    # compute VX and VL -
    kEX = max_transcription_rate*(3600/characteristic_gene_length)  # hr^-1
    kEL = max_translation_rate*(3600/characteristic_prot_length)    # hr^-1
    VX = kEX*RX_concentration                                       # μM/hr
    VL = kEL*RL_concentration                                       # μM/hr

    # compute saturation coefficient for X and L -
    m = parse(Float64,simulation_json_tree["biophysical_parameters"]["elongation_slope"]["value"]) # μM
    KX = m                  # μM
    KL = 190*KX             # μM

    # get the induction parameters -
    induction_parameter_dictionary = Dict{Symbol,Any}()
    K1 = parse(Float64,simulation_json_tree["simulation_paramaters"]["induction_parameters"]["K1"])
    K2 = parse(Float64,simulation_json_tree["simulation_paramaters"]["induction_parameters"]["K2"])
    n = parse(Float64,simulation_json_tree["simulation_paramaters"]["induction_parameters"]["coop_parameter"])
    k_binding = parse(Float64,simulation_json_tree["simulation_paramaters"]["induction_parameters"]["k_binding"])
    induction_parameter_dictionary[:K1] = K1
    induction_parameter_dictionary[:K2] = K2
    induction_parameter_dictionary[:coop_parameter] = n
    induction_parameter_dictionary[:k_binding] = k_binding
    induction_parameter_dictionary[:inducer_concentration] = 0.0 # default value

    # == DO NOT EDIT BELOW THIS LINE ========================================== #
    # parameter_dictionary[:stoichiometric_matrix] = balanced_stoichiometrix_matrix
    # parameter_dictionary[:flux_bounds_array] = flux_bounds_array
    # parameter_dictionary[:species_bounds_array] = species_bounds_array
    # parameter_dictionary[:objective_coefficient_array] = objective_coefficient_array

    # stuff for bounds -
    parameter_dictionary[:gene_concentration] = gene_concentration
    parameter_dictionary[:length_factor_transcription] = length_factor_transcription
    parameter_dictionary[:length_factor_translation] = length_factor_translation
    parameter_dictionary[:RX_concentration] = RX_concentration
    parameter_dictionary[:RL_concentration] = RL_concentration
    parameter_dictionary[:X_tau_factor] = X_tau_factor
    parameter_dictionary[:L_tau_factor] = L_tau_factor
    parameter_dictionary[:kdX] = kdX
    parameter_dictionary[:kdL] = kdL
    parameter_dictionary[:KX] = KX
    parameter_dictionary[:KL] = KL
    parameter_dictionary[:kdL] = kdL
    parameter_dictionary[:VX] = VX
    parameter_dictionary[:VL] = VL
    parameter_dictionary[:induction_parameter_dictionary] = induction_parameter_dictionary

    return parameter_dictionary
    # ========================================================================= #
end
# -- PRIVATE HELPER METHODS ------------------------------------------------------------------------------- #
