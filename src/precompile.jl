using PrecompileTools: @compile_workload, @setup_workload

function _precompile_workload()
    sbml = """
    <?xml version="1.0" encoding="UTF-8"?>
    <sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">
      <model id="precompile_model">
        <listOfCompartments>
          <compartment id="cell" size="1" constant="true"/>
        </listOfCompartments>
        <listOfSpecies>
          <species id="x" compartment="cell" initialAmount="1"
                   hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
        </listOfSpecies>
        <listOfParameters>
          <parameter id="k" value="0.5" constant="true"/>
        </listOfParameters>
        <listOfReactions>
          <reaction id="decay" reversible="false">
            <listOfReactants>
              <speciesReference species="x" stoichiometry="1" constant="true"/>
            </listOfReactants>
            <kineticLaw>
              <math xmlns="http://www.w3.org/1998/Math/MathML">
                <apply><times/><ci>k</ci><ci>x</ci></apply>
              </math>
            </kineticLaw>
          </reaction>
        </listOfReactions>
      </model>
    </sbml>
    """

    checksupport_string(sbml)
    model = readSBMLFromString(sbml)
    ReactionSystem(model)
    ODESystem(model)
    return nothing
end

@setup_workload begin
    @compile_workload begin
        _precompile_workload()
    end
end
