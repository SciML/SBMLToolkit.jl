using SBML

sbml_fn = joinpath(@__DIR__, "data/simpleModel.sbml")

m_ = readSBML(sbml_fn, doc -> begin
    set_level_and_version(3, 2)(doc)
    convert_simplify_math(doc)
end)

m = deepcopy(m_)
m2 = SBML.fix_simbio_names(m)

@test keys(m2.species) == Set(["unnamed₊A", "unnamed₊B"])
