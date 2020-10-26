"""
    File name: dfoundations_test.py
    Author: Rob Swart (Mobilis TBI)
    Additional Notes: Niels van der Zee (RHDHV) --> See #NZE
    Date last modified: 14/09/2020
    Script to test D-Foundations Geolib API
    Construction type: prefab and tubex bearing piles
    Calculation type: Bearing piles EC7-NL, complete verification calculation
"""
import time
start = time.time()
import geolib as gl
from geolib import DFoundationsModel
from geolib.geometry import Point
from geolib.models.dfoundations import profiles
from geolib.models.dfoundations.dfoundations_model import BearingPilesModel, CalculationOptions
from geolib.soils import Soil
from geolib.models.dfoundations import piles
import pandas as pd
from pathlib import Path
from pydantic.color import Color

def init_Model(pile_tip_level, Fd, Fk):

    df = DFoundationsModel()
    # Set model: should be done at beginning of workflow because all soils will be replaced by defaults! So custom soils will get lost
    # Bearing piles NL, complete verification calculation
    # todo Test other calculation types and tension piles
    model_options = BearingPilesModel(
        is_rigid=False, factor_xi3=1.2, factor_xi4=1.2
    )
    calculation_options = CalculationOptions(
        calculationtype=gl.models.dfoundations.dfoundations_model.CalculationType.VERIFICATION_COMPLETE,
        cpt_test_level=pile_tip_level,
    )
    df.set_model(model_options, calculation_options)
    return df

def add_soil(df, soil_name,ysat,yunsat,phi,su,d50,soil_type):
    print(soil_name,ysat,yunsat,phi,su,d50,soil_type)
    new_soil = Soil()
    new_soil.name = str(soil_name)
    new_soil.soil_weight_parameters.saturated_weight = float(ysat)
    new_soil.soil_weight_parameters.unsaturated_weight = float(yunsat)
    new_soil.mohr_coulomb_parameters.friction_angle = float(phi)
    new_soil.undrained_parameters.undrained_shear_strength = float(su)  # Required, but I don't understand why
    new_soil.soil_classification_parameters.d_50 = float(d50)
    new_soil.soil_type_nl = soil_type  # Options should be explained in manual
    # sand1.color= int(10871211) #color to be added later
    df.add_soil(new_soil)

def add_cpt(cptname, GL, data, timeorder):
   # Define CPTs
    cpt = profiles.CPT(
        cptname=cptname,
        groundlevel=GL,
        measured_data=data,
        timeorder_type=timeorder,)
    return cpt


def add_profile(df,cptname, GL, data, timeorder, layers, excavation, x, y , phreatic_level, pile_tip_level, OCR, ToPSF, BoNSF, Settlement):

    cpt = add_cpt(cptname, GL, data, timeorder)
    profile1 = profiles.Profile(
            name="CPT1",
            cpt=cpt, #this refers to cpt created in add_cpt (check how to bind these later
            location=profiles.Point(x=x, y=y),
            layers=layers ,
            phreatic_level= phreatic_level,
            pile_tip_level= pile_tip_level,
            overconsolidation_ratio= OCR,
            top_of_positive_skin_friction=ToPSF ,
            bottom_of_negative_skin_friction=BoNSF ,
            expected_ground_level_settlement= Settlement,
            excavation=excavation,
        )
    df.add_profile(profile1)

def new_pile_type(Piles):

    new_parent_pile = dict(
        pile_name=Piles["pile_name"],
        pile_type=Piles["pile_type"],
        execution_factor_sand_gravel=Piles["execution_factor_sand_gravel"],
        # Only first option (with STANDARD) seems to work (in this case, no factor input necessary)
        pile_type_for_execution_factor_clay_loam_peat=Piles["pile_type_for_execution_factor_clay_loam_peat"],
        # Better name: preset_pile_class_factor_shaft_clay_loam_peat

        # pile_type_for_execution_factor_clay_loam_peat=piles.BasePileTypeForClayLoamPeat.USER_DEFINED,     # Better name: preset_pile_class_factor_shaft_clay_loam_peat
        # execution_factor_clay_loam_peat=0.01, # Better name: pile_class_factor_shaft_clay_loam_peat -> doesn't seem to work

        pile_class_factor=Piles["pile_class_factor"],  # Better name: pile_class_factor_tip
        load_settlement_curve = Piles["load_settlement_curve"],
        user_defined_pile_type_as_prefab= Piles["user_defined_pile_type_as_prefab"],
        use_manual_reduction_for_qc= Piles["use_manual_reduction_for_qc"],
        elasticity_modulus= Piles["elasticity_modulus"],
        characteristic_adhesion= Piles["characteristic_adhesion"],
        overrule_pile_tip_shape_factor= Piles["overrule_pile_tip_shape_factor"],
        overrule_pile_tip_cross_section_factors= Piles["overrule_pile_tip_cross_section_factors"],
    )
    new_pile = piles.BearingRectangularPile(**new_parent_pile, **Piles["geometry"])

    return new_pile


def add_piles(df, Pileplan, pile_type):
    # Define pile properties
    Piles = {}
    #Piles = [piletype, X,Y, Pile_head, surcharge, limit_state_str, limit_state_service]
    #Fill out later, maybe a dict???
    counter = 0
    for i in range(0,len(Pileplan)):
        location = piles.BearingPileLocation(
            pile_name=counter+1, #perhaps number through? with some counter?
            point=Point(x=Pileplan.loc[i]["X"], y=Pileplan.loc[i]["Y"]),
            pile_head_level=Pileplan.loc[i]["Pile_head_level"],
            surcharge=Pileplan.loc[i]["Surcharge"],
            limit_state_str=Pileplan.loc[i]["limit_state_str"],
            limit_state_service=Pileplan.loc[i]["limit_state_service"],
        )
        df.add_pile_if_unique(pile_type, location)
        counter += 1

if __name__ == '__main__':



    Soils = pd.DataFrame(columns=["soil_name","ysat","yunsat","phi","su","d50","soil_type"])
    Soils.loc[0] = ["Sand type 1", 19, 17, 30, 1, 0.2, 1]
    Soils.loc[1]= ["Sand type 2", 21, 19, 40, 1, 0.2, 1]
    Soils.loc[2] = ["Clay type 1", 15, 15, 20, 1, 0.2, 3]


    Calculation = init_Model(-15, 800, 500) #Depends on the calculation type as well.

    for i in range(0,len(Soils)):
        add_soil(Calculation, Soils.loc[i]["soil_name"],  Soils.loc[i]["ysat"],  Soils.loc[i]["yunsat"],
                 Soils.loc[i]["phi"],  Soils.loc[i]["su"], Soils.loc[i]["d50"], Soils.loc[i]["soil_type"])

    cptname = "cpt1"
    cpt_data = [
            {"z": 0.0, "qc": 0.1},
            {"z": -0.10, "qc": 0.5},
            {"z": -0.20, "qc": 2.0},
            {"z": -0.30, "qc": 3.0},
            {"z": -0.40, "qc": 5.0},
            {"z": -10, "qc": 1.0},
            {"z": -15, "qc": 5.0},
            {"z": -25, "qc": 13.0},
            {"z": -30, "qc": 22.0},
        ]

    cpt_ground_level = 0.5
    cpt_timeorder = profiles.TimeOrderType.CPT_EXCAVATION_INSTALL

    layers = [
        {
            "material": "Sand type 1",
            "top_level": 0.0,
            "excess_pore_pressure_top": 0.0,
            "excess_pore_pressure_bottom": 0.0,
            "ocr_value": 1.0,
            "reduction_core_resistance": 0,
        },
        {
            "material": "Clay type 1",
            "top_level": -3,
            "excess_pore_pressure_top": 0.0,
            "excess_pore_pressure_bottom": 0.0,
            "ocr_value": 1.0,
            "reduction_core_resistance": 0,
        },
        {
            "material": "Sand type 2",
            "top_level": -5,
            "excess_pore_pressure_top": 0.0,
            "excess_pore_pressure_bottom": 0.0,
            "ocr_value": 1.0,
            "reduction_core_resistance": 0,
        },
    ]
    #There seems to be no distinction between filling out the distance edge to excavation or not (How to turn on begemann?)
    excavation = profiles.Excavation(excavation_level=-5,
                                     distance_edge_pile_to_excavation_boundary = 2.00)
    x = 0
    y = 5.0
    phreatic_level = -0.5
    pile_tip_level = -15
    OCR =  1
    ToPSF = -6
    BoNSF = 0
    Settlement = 0

    add_profile(Calculation, cptname, cpt_ground_level, cpt_data, cpt_timeorder,
                layers, excavation, x, y, phreatic_level, pile_tip_level, OCR, ToPSF,
                BoNSF, Settlement)

    Piles = pd.DataFrame(columns=["geometry", #Create a dict (see documentation for more information
                                  "pile_name", #User defined name
                                  "pile_type", #See documentation for possibilities
                                  "execution_factor_sand_gravel", #either user defined or choose piletype (see documentation)
                                  "pile_type_for_execution_factor_clay_loam_peat", #either user defined or choose piletype (see documentation)
                                  "execution_factor_clay_loam_peat", #either user defined or choose piletype (see documentation)
                                  "pile_class_factor", #either user defined or choose piletype (see documentation)
                                  "load_settlement_curve", #see documentation for options
                                  "user_defined_pile_type_as_prefab", #either user defined or choose piletype (see documentation)
                                  "use_manual_reduction_for_qc", #Boolean
                                  "elasticity_modulus", #value
                                  "characteristic_adhesion", #value
                                  "overrule_pile_tip_shape_factor",#Boolean
                                  "overrule_pile_tip_cross_section_factors"]) #Boolean


    Piles.loc[0] = [dict(base_width=0.45, base_length=0.45),"Prefab",piles.BasePileType.USER_DEFINED_VIBRATING,0.01,
                piles.BasePileTypeForClayLoamPeat.STANDARD,0.01, 0.7,piles.LoadSettlementCurve.ONE,False,False,
                2*10**7,0,False,False]


    Prefab = new_pile_type(Piles.loc[0])
    Pileplan = pd.DataFrame(columns=["X","Y","Pile_head_level","Surcharge","limit_state_str","limit_state_service"])
    for i in range(0,9):
        Pileplan.loc[i] = [i*10,i*10,0,0,800,500]

    add_piles(Calculation, Pileplan,Prefab)

    print('serializing')
    dfoundations_file = Path("dfoundations_test_1.foi")
    Calculation.serialize(dfoundations_file)
    Calculation.execute()

    output_dict = Calculation.output.dict()



    Rd = output_dict['verification_results']['nen_pile_results']['max_bearing_capacity_foundation']
    Fd = 800
    Fk = 500
    print("### Results for pile tip at " + str(pile_tip_level) + " m NAP ###" +
          "\n\n## Design bearing capacity ##" +
          "\nFd= " + str(Fd) + " kN" +
          "\nRd= " + str(round(Rd)) + " kN" +
          "\nuc= " + str(round(Fd / Rd, 2)))

    s_neg1b = output_dict['verification_results']['nen_pile_results']['sneg1b']
    s_b1b = output_dict['verification_results']['nen_pile_results']['spunt_d_1b']
    s_el_d1b = output_dict['verification_results']['nen_pile_results']['sel_d1b']
    s_21b = output_dict['verification_results']['nen_pile_results']['s2_d1b']
    s_d1b = s_neg1b + s_b1b + s_el_d1b + s_21b

    print("\n## Settlement STR/GEO ##" +
          "\ns_neg= " + str(s_neg1b) + " m" +
          "\ns_b= " + str(s_b1b) + " m" +
          "\ns_el;d= " + str(s_el_d1b) + " m" +
          "\ns2= " + str(s_21b) + " m" +
          "\n----------------" +
          "\ns_d_tot= " + str(round(s_d1b, 5)) + " m")

    s_neg2 = output_dict['verification_results']['nen_pile_results']['sneg2']
    s_b2 = output_dict['verification_results']['nen_pile_results']['spunt_d_2']
    s_el_d2 = output_dict['verification_results']['nen_pile_results']['sel_d2']
    s_22 = output_dict['verification_results']['nen_pile_results']['s2_d2']
    s_d2 = s_neg2 + s_b2 + s_el_d2 + s_22

    print("\n## Settlement SLS ##" +
          "\ns_neg= " + str(s_neg2) + " m" +
          "\ns_b= " + str(s_b2) + " m" +
          "\ns_el;d= " + str(s_el_d2) + " m" +
          "\ns2= " + str(s_22) + " m" +
          "\n----------------" +
          "\ns_k_tot= " + str(round(s_d2, 5)) + " m")

    print("\n## Axial pile stiffness (characteristic) ##" +
          "\nFk/s_k_tot= " + str(Fk) + "/" + str(round(s_d2, 5)) + "= " + str(round(Fk / s_d2 / 1000)) + " MN/m")

    print("\n## Max shaft and point resistance ##")
    max_shaft_and_point = output_dict['verification_results']['nen_pile_results']['max_shaft_and_point']['data']
    print(pd.DataFrame(max_shaft_and_point).reindex(columns=['CptIndex', 'Grenstoestand', 'MaxShaft', 'MaxPoint']))

    print("\n## Used CPTs ##")
    used_CPTs = output_dict['verification_results']['nen_pile_results']['cpts']['data']
    print(pd.DataFrame(used_CPTs).reindex(
        columns=['index', 'PPN', 'HNK', 'HPK', 'XCoordinate', 'YCoordinate', 'CptName']))
