module PseudoPotentials
using AtomicLevels
import AtomicLevels: spectroscopic_label, AbstractOrbital
using AtomicPotentials
using PrettyTables

struct PseudoPotential{T,relativistic} <: AbstractPotential{T}
    name::String
    gst_config::Configuration{Orbital{Int}}
    Q::Int
    Vℓ::Vector{Matrix{T}} # ℓ ∈ 0:ℓmax
    Vℓ′::Vector{Matrix{T}} # ℓ′ ∈ 1:ℓmax′
    reference::String
end
AtomicPotentials.charge(pp::PseudoPotential) = num_electrons(pp.gst_config)
AtomicPotentials.ground_state(pp::PseudoPotential) = pp.gst_config

function Base.show(io::IO, pp::PseudoPotential{T,relativistic}) where {T,relativistic}
    if relativistic
        write(io, "Relativistic pseudo-potential")
    else
        write(io, "Pseudo-potential")
    end
    write(io, " for $(pp.name) ($(pp.gst_config)), Z = $(charge(pp))")
end

function Base.show(io::IO, ::MIME"text/plain", pp::PseudoPotential)
    show(io, pp)
    println(io)
    ℓmax = length(pp.Vℓ)-1
    ℓmax′ = length(pp.Vℓ′)
    print(io, "Long-range Q = $(pp.Q), ℓ ∈ 0:$(ℓmax)")
    ℓmax′ > 0 && print(io, ", ℓ′ ∈ 1:$(length(pp.Vℓ′))")
    println(io, "\nData from \"$(pp.reference)\"")
    ℓk = map(enumerate(0:ℓmax)) do (i,ℓ)
        nk = size(pp.Vℓ[i],1)
        [[spectroscopic_label(ℓ);repeat([""],nk-1)] 1:nk]
    end |> v -> vcat(v...)
    ℓ′k = map(enumerate(1:ℓmax′)) do (i,ℓ′)
        nk = size(pp.Vℓ′[i],1)
        [[spectroscopic_label(ℓ′);repeat([""],nk-1)] 1:nk]
    end |> v -> vcat(v...)
    data = [ℓk vcat(pp.Vℓ...)]
    headers = ["ℓ", "k", "n", "β", "B"]
    if ℓmax′ > 0
        data = [data vcat(repeat([""], size(pp.Vℓ[1],1), 5),
                          [ℓ′k vcat(pp.Vℓ′...)],
                          repeat([""], 1, 5))]
        headers = vcat(headers, ["ℓ′", "k", "n", "β", "B"])
    end
    pretty_table(io, data, headers)
end

#=
\[V_{\mathrm{pp}}(\vec{r})=
-\frac{Q}{r}+
\sum_{\ell jk}
B_{\ell j k}
r^{n_k-2}
\exp(-\beta_{\ell j k} r^2)
\mathcal{P}_{\ell j}\]
=#

function (pp::PseudoPotential{T,true})(orb::RelativisticOrbital, r::AbstractVector{T}) where T
    V = -pp.Q./r
    κ = orb.κ
    if isempty(pp.Vℓ′) && κ > 0
        κ = AtomicLevels.flip_j(orb).κ
    end
    data = if κ < 0
        -κ > length(pp.Vℓ) && return V
        pp.Vℓ[-κ]
    else
        κ > length(pp.Vℓ′) && return V
        pp.Vℓ′[κ]
    end
    for k in 1:size(data,1)
        (n,β,B) = data[k,:]
        V += B * r.^(n-2).*exp.(-β*r.^2)
    end
    V
end

function (pp::PseudoPotential{T,true})(orb::SpinOrbital, r::AbstractVector{T}) where T
    V = zeros(T, length(r))
    for ((rorb,mj),coeff) in jj2lsj(orb)
        V += coeff*pp(rorb, r)
    end
    V
end

function (pp::PseudoPotential{T,false})(orb::SpinOrbital, r::AbstractVector{T}) where T
    V = -pp.Q./r
    ℓ = orb.orb.ℓ
    ℓ >= length(pp.Vℓ) && return V
    data = pp.Vℓ[ℓ+1]
    for k in 1:size(data,1)
        (n,β,B) = data[k,:]
        V += B * r.^(n-2).*exp.(-β*r.^2)
    end

    V
end

(pp::PseudoPotential{T})(orb::AbstractOrbital, r::T) where T = pp(orb, [r])[1]

function parse_pseudopotential_line(::Type{T}, pp_line::AbstractString) where T
    data = split(pp_line, ";")
    num_terms = parse(Int, data[1])
    length(data) == num_terms+2 || throw(ArgumentError("Unrecognizable potential string $(pp_line)"))
    map(1:num_terms) do i
        term = split(strip(data[i+1]),",")
        parse.(Ref(T), term)'
    end |> p -> vcat(p...)
end

function parse_pseudopotential(::Type{T}, pp_string::String, relativistic::Bool) where T
    lines = split(pp_string, "\n")
    gst_config = parse(Configuration{Orbital},
                       lstrip(lines[1], ['!', ' ']))
    i = findfirst(l -> first(l) != '!', lines)
    i === nothing && error("Could not find start of data")
    header = split(split(lines[i], ";")[1], ",")
    kind = header[1]

    element = String(header[2])
    Z = element_number(Symbol(element))
    nelec = num_electrons(gst_config)
    Z == nelec ||
        throw(ArgumentError("Atom number Z = $Z does not match number of electrons of ground state configuration $(gst_config) => $(nelec) electrons"))

    Q̃ = parse(Int, header[3])
    ncore = num_electrons(core(gst_config))
    Q̃ == ncore ||
        throw(ArgumentError("Charge modelled by pseudopotential, Q̃ = $(Q̃), does not match number of core electrons in ground state configuration: $(core(gst_config)) => $(ncore) electrons"))

    Q = Z-Q̃ # Long-range charge of potential

    ℓmax = parse(Int, header[4])
    ℓmax′ = parse(Int, header[5])
    ii = i
    Vℓ = map(vcat(ℓmax,0:ℓmax-1)) do ℓ
        ℓ => parse_pseudopotential_line(T, lines[ii+=1])
    end |> v -> last.(sort(v, by=first))
    Vℓ′ = map(1:ℓmax′) do ℓ′
        parse_pseudopotential_line(T, lines[ii+=1])
    end |> Vector{Matrix{T}}
    reference = map(lines[ii+2:end]) do line
        lstrip(line, ['!', '[', ']', ' ', '0':'9'...])
    end |> l -> join(l, "\n")
    PseudoPotential{T,relativistic}(element, gst_config, Q, Vℓ, Vℓ′, reference)
end

macro PP_str(pp_string)
    parse_pseudopotential(Float64, pp_string, false)
end

macro RPP_str(pp_string)
    parse_pseudopotential(Float64, pp_string, true)
end

# * Misc pseudopotentials
# ** Neon

NeonHF = PP"""! [He] 2s2 2p6
!  Q=8., MEFIT, HF, Ref 22.
ECP,Ne,2,3,0;
1; 2,1.000000,0.00000000;
2; 2,32.029855,112.52822460; 2,12.294136,28.44589604;
2; 2,21.525406,-10.93982872; 2,13.091176,3.39426049;
1; 2,0.850385,-0.16409678;
! References:
! [22] A. Nicklass, M. Dolg, H. Stoll, H. Preuss, J. Chem. Phys. 102, 8942 (1995)."""

NeonWB = PP"""! [He] 2s2 2p6
!  Q=8., MEFIT, WB, Ref 22.
ECP,Ne,2,3,0;
1; 2,1.000000,0.00000000;
2; 2,31.860162,112.52543566; 2,12.362219,28.30083454;
2; 2,21.508034,-11.12658543; 2,12.910447,3.38754919;
1; 2,0.850385,-0.18408921;
! References:
! [22] A. Nicklass, M. Dolg, H. Stoll, H. Preuss, J. Chem. Phys. 102, 8942 (1995)."""

# ** Argon

ArgonHF = PP"""! [Ne] 3s2 3p6
!  Q=8., MEFIT, HF, Ref 22.
ECP,Ar,10,4,0;
1; 2,1.000000,0.00000000;
2; 2,10.291970,68.69327818; 2,3.947937,24.48814203;
2; 2,5.388907,27.59672144; 2,2.704463,4.13654627;
2; 2,8.084705,-8.14025335; 2,4.018469,-1.66241708;
1; 2,5.295301,-3.40922653;
! References:
! [22] A. Nicklass, M. Dolg, H. Stoll, H. Preuss, J. Chem. Phys. 102, 8942 (1995)."""

ArgonWB = PP"""! [Ne] 3s2 3p6
!  Q=8., MEFIT, WB, Ref 22.
ECP,Ar,10,4,0;
1; 2,1.000000,0.00000000;
2; 2,10.261721,68.66778801; 2,3.952725,24.04276629;
2; 2,5.392714,27.73076331; 2,2.699967,4.04545904;
2; 2,8.086235,-8.13747696; 2,4.016632,-1.66452808;
1; 2,5.208459,-3.40009845;
! References:
! [22] A. Nicklass, M. Dolg, H. Stoll, H. Preuss, J. Chem. Phys. 102, 8942 (1995)."""

# ** Krypton

KryptonDF = RPP"""! [Ne] 3s2 3p6 3d10 4s2 4p6
!  Q=26., MEFIT, MCDHF+Breit, Ref 36.
ECP,Kr,10,4,3;
1; 2,1.000000,0.000000;
3; 2,70.034410,49.953835; 2,31.895971,369.978238; 2,7.353728,10.054544;
4; 2,47.265183,99.101833; 2,46.664268,198.232524; 2,23.008133,28.204670; 2,22.100979,56.516792;
6; 2,50.768165,-18.588877; 2,50.782228,-27.875398; 2,15.479497,-0.293411; 2,15.559383,-0.469399; 2,2.877154,0.068881; 2,1.991196,0.132270;
2; 2,15.437567,-1.223389; 2,22.055742,-2.991233;
4; 2,47.265183,-198.203666; 2,46.664268,198.232524; 2,23.008133,-56.409339; 2,22.100979,56.516792;
6; 2,50.768165,18.588877; 2,50.782228,-18.583599; 2,15.479497,0.293411; 2,15.559383,-0.312933; 2,2.877154,-0.068881; 2,1.991196,0.088180;
2; 2,15.437567,0.815593; 2,22.055742,-1.495617;
! References:
! [36] K.A. Peterson, D. Figgen, E. Goll, H. Stoll, M. Dolg, J. Chem. Phys. 119, 11113 (2003).
"""

# ** Xenon

XenonHF = PP"""! [Ar] 3d10c 4s2c 4p6c 4d10c 5s2 5p6
!  Q=8., MEFIT, HF, Ref 22.
ECP,Xe,46,4,0;
1; 2,1.000000,0.00000000;
2; 2,3.815600,122.76231371; 2,1.878604,8.20135456;
2; 2,3.003078,68.75731963; 2,1.283819,3.64084871;
2; 2,2.027610,23.08319830; 2,0.830435,3.17298823;
2; 2,6.246157,-48.30201724; 2,1.564998,-6.91572892;
! References:
! [22] A. Nicklass, M. Dolg, H. Stoll, H. Preuss, J. Chem. Phys. 102, 8942 (1995)."""

XenonWB = PP"""! [Ar] 3d10c 4s2c 4p6c 4d10c 5s2 5p6
!  Q=8., MEFIT, WB, Ref 22; CPP: alpha=0.8345;delta=0.97391;ncut=1.
ECP,Xe,46,5,0;
1; 2,1.000000,0.00000000;
2; 2,3.940263,122.76382934; 2,2.277264,8.30885115;
2; 2,3.028373,68.82300437; 2,1.394319,3.64674223;
2; 2,2.122605,23.65207854; 2,0.798669,3.25844113;
2; 2,6.164360,-47.70319876; 2,1.542374,-6.54113991;
1; 2,1.847892,-7.10585060;
! References:
! [22] A. Nicklass, M. Dolg, H. Stoll, H. Preuss, J. Chem. Phys. 102, 8942 (1995)."""

XenonDF = RPP"""! [Ar] 3d10c 4s2 4p6 4d10 5s2 5p6
!  Q=26., MEFIT, MCDHF+Breit, Ref 36.
ECP,Xe,28,4,3;
1; 2,1.000000,0.000000;
3; 2,40.005184,49.997962; 2,17.812214,281.013303; 2,9.304150,61.538255;
4; 2,15.701772,67.439142; 2,15.258608,134.874711; 2,9.292184,14.663300; 2,8.559003,29.354730;
6; 2,15.185600,35.436908; 2,14.284500,53.195772; 2,7.121889,9.046232; 2,6.991963,13.223681; 2,0.623946,0.084853; 2,0.647284,0.044155;
4; 2,20.881557,-23.089295; 2,20.783443,-30.074475; 2,5.253389,-0.288227; 2,5.361188,-0.386924;
4; 2,15.701772,-134.878283; 2,15.258608,134.874711; 2,9.292184,-29.326600; 2,8.559003,29.354730;
6; 2,15.185600,-35.436908; 2,14.284500,35.463848; 2,7.121889,-9.046232; 2,6.991963,8.815787; 2,0.623946,-0.084853; 2,0.647284,0.029436;
4; 2,20.881557,15.392863; 2,20.783443,-15.037237; 2,5.253389,0.192151; 2,5.361188,-0.193462;
! References:
! [36] K.A. Peterson, D. Figgen, E. Goll, H. Stoll, M. Dolg, J. Chem. Phys. 119, 11113 (2003)."""

# ** Radon

RadonDF = RPP"""! [Kr] 4d10c 4f14c 5s2 5p6 5d10 6s2 6p6
!  Q=26., MEFIT, MCDHF+Breit, Ref 36.
ECP,Rn,60,5,4;
1; 2,1.000000,0.000000;
3; 2,30.151242,49.965551; 2,14.521241,283.070000; 2,8.052038,62.002870;
4; 2,11.009942,71.969119; 2,9.617625,143.860559; 2,7.336008,4.714761; 2,6.406253,9.013065;
4; 2,8.369220,36.368365; 2,8.116975,54.551761; 2,5.353656,9.634487; 2,5.097212,14.387902;
4; 2,6.348571,21.797290; 2,6.295949,28.946805; 2,2.882118,1.447365; 2,2.908048,2.177964;
4; 2,11.015205,-25.228095; 2,10.909853,-30.566121; 2,3.482761,-0.574800; 2,3.600418,-0.779538;
4; 2,11.009942,-143.938238; 2,9.617625,143.860559; 2,7.336008,-9.429523; 2,6.406253,9.013065;
4; 2,8.369220,-36.368365; 2,8.116975,36.367841; 2,5.353656,-9.634487; 2,5.097212,9.591935;
4; 2,6.348571,-14.531526; 2,6.295949,14.473402; 2,2.882118,-0.964910; 2,2.908048,1.088982;
4; 2,11.015205,12.614048; 2,10.909853,-12.226448; 2,3.482761,0.287400; 2,3.600418,-0.311815;
! References:
! [36] K.A. Peterson, D. Figgen, E. Goll, H. Stoll, M. Dolg, J. Chem. Phys. 119, 11113 (2003)."""

export @PP_str, @pc_str, AbstractPotential, charge, ground_state

end # module
