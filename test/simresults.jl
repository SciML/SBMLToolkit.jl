using SBMLToolkitTestSuite
using Test

const case_ids = [7,  # boundary_condition
    22,  # non-integer stoichiometry
    23,  # species with constant=boundaryCondition="true"
    41,  # events
    140,  # compartment size overridden with assignmentRule
    170,  # Model using parameters and rules only
    325,  # One reactions and two rate rules with four species in a 2D compartment
    679  # Initial value calculated by assignmentRule in compartment of non-unit size    # 1208, # Non-hOSU species with rateRule in variable compartment -> require MTK fix.
]

const logdir = joinpath(@__DIR__, "logs")
ispath(logdir) && rm(logdir, recursive = true)
mkdir(logdir)

df = verify_all(case_ids, logdir)

for i in 1:length(case_ids)
    @test sum(Vector(df[i, ["expected_err", "res"]])) == 1
end
