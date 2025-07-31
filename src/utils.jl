# Conversion to symbolics
const IV = default_t()
const D = default_time_deriv()
symbolicsRateOf(x) = D(x)
const symbolics_mapping = Dict(SBML.default_function_mapping...,
    "rateOf" => symbolicsRateOf)

function interpret_as_num(x::SBML.Math, model::SBML.Model)
    SBML.interpret_math(x;
        map_apply = (x::SBML.MathApply,
            interpret::Function) -> Num(symbolics_mapping[x.fn](interpret.(x.args)...)),
        map_const = (x::SBML.MathConst) -> Num(SBML.default_constants[x.id]),
        map_ident = x -> map_symbolics_ident(x, model),
        map_lambda = (_,
            _) -> throw(ErrorException("Symbolics.jl does not support lambda functions")),
        map_time = (x::SBML.MathTime) -> IV,
        map_value = (x::SBML.MathVal) -> x.val,
        map_avogadro = (x::SBML.MathAvogadro) -> SBML.default_constants["avogadro"])
end

"""
Get dictionary to change types in kineticLaw
"""
function get_substitutions(model)
    u0map, parammap = get_mappings(model)
    subsdict = Dict()
    for item in first.(u0map)
        k = create_var(string(item.f.name))
        subsdict[k] = item
    end
    for item in first.(parammap)
        k = create_var(string(item.name))
        subsdict[k] = item
    end
    subsdict
end

function map_symbolics_ident(x::SBML.Math, model::SBML.Model)
    k = x.id
    category = k in keys(model.species) ? :species :
               k in keys(model.parameters) ? :parameter :
               k in keys(model.compartments) ? :compartment :
               error("Unknown category for $k")
    if k in keys(model.species)
        v = model.species[k]
        if v.constant == true
            var = create_param(k; isconstantspecies = true)
        else
            var = create_var(k, IV;
                isbcspecies = has_rule_type(k, model, SBML.RateRule) ||
                              has_rule_type(k, model, SBML.AssignmentRule) ||
                              (has_rule_type(k, model, SBML.AlgebraicRule) &&
                               (all([netstoich(k, r) == 0
                                     for r in values(model.reactions)]) ||
                                v.boundary_condition == true)))  # To remove species that are otherwise defined
        end
    elseif k in keys(model.parameters)
        v = model.parameters[k]
        if v.constant == false &&
           (SBML.seemsdefined(k, model) || is_event_assignment(k, model))
            var = create_var(k, IV; isbcspecies = true)
        elseif v.constant == true && isnothing(v.value)  # Todo: maybe add this branch also to model.compartments
            var = create_param(k)
        else
            var = create_param(k)
        end
    elseif k in keys(model.compartments)
        v = model.compartments[k]
        if v.constant == false && SBML.seemsdefined(k, model)
            var = create_var(k, IV; isbcspecies = true)
        else
            var = create_param(k)
        end
    else
        error("$k must be in the model species, parameters, or compartments.")
    end
    Num(var)
end

function create_var(x; isbcspecies = false)
    sym = Symbol(x)
    Symbolics.unwrap(first(@species $sym [isbcspecies = isbcspecies]))
end
function create_var(x, iv; isbcspecies = false, irreducible = false)
    sym = Symbol(x)
    Symbolics.unwrap(first(@species $sym(iv) [
        isbcspecies = isbcspecies,
        irreducible = irreducible
    ]))
end
function create_param(x; isconstantspecies = false)
    sym = Symbol(x)
    Symbolics.unwrap(first(@parameters $sym [isconstantspecies = isconstantspecies]))
end

function has_rule_type(id::String, m::SBML.Model, T::Type{<:SBML.Rule})
    T == SBML.AlgebraicRule &&
        return any(SBML.isfreein(id, r.math) for r in m.rules if r isa SBML.AlgebraicRule)
    any(r.variable == id for r in m.rules if r isa T)
end

const importdefaults = doc -> begin
    set_level_and_version(3, 2)(doc)
    convert_promotelocals_expandfuns(doc)
end

function create_symbol(k::String, model::SBML.Model)
    if k in keys(model.species)
        v = model.species[k]
        if v.constant == true
            sym = create_param(k; isconstantspecies = true)
        else
            sym = create_var(k, IV;
                isbcspecies = has_rule_type(k, model, SBML.RateRule) ||
                              has_rule_type(k, model, SBML.AssignmentRule) ||
                              (has_rule_type(k, model, SBML.AlgebraicRule) &&
                               (all([netstoich(k, r) == 0
                                     for r in values(model.reactions)]) ||
                                v.boundary_condition == true)))  # To remove species that are otherwise defined
        end
    elseif k in keys(model.parameters)
        v = model.parameters[k]
        if v.constant == false &&
           (SBML.seemsdefined(k, model) || is_event_assignment(k, model))
            sym = create_var(k, IV; isbcspecies = true)
        else
            sym = create_param(k)
        end
    elseif k in keys(model.compartments)
        v = model.compartments[k]
        if v.constant == false && SBML.seemsdefined(k, model)
            sym = create_var(k, IV; isbcspecies = true)
        else
            sym = create_param(k)
        end
    end
    sym
end
