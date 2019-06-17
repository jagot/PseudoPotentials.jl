module PseudoPotentials
using AtomicLevels
import AtomicLevels: spectroscopic_label, AbstractOrbital
using AtomicPotentials
using AngularMomentumAlgebra
import AngularMomentumAlgebra: jmⱼ
using WignerSymbols
using PrettyTables

abstract type AbstractPseudoPotential{T} <: AbstractPotential{T} end

include("gaussian_expansions.jl")
include("nonrelativistic_pseudopotentials.jl")
include("relativistic_pseudopotentials.jl")
include("parse.jl")
include("misc_pseudopotentials.jl")

export @PP_str, @pc_str,
    AbstractPotential, PseudoPotential, RelativisticPseudoPotential,
    charge, ground_state, spin_orbit_potential

end # module
