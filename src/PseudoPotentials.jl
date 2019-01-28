module PseudoPotentials
using AtomicLevels
import AtomicLevels: spectroscopic_label
using AtomicPotentials
using PrettyTables

struct PseudoPotential{T} <: AbstractPotential{T}
    name::String
    gst_config::Configuration{Orbital{Int}}
    Q::Int
    Vℓ::Vector{Matrix{T}} # ℓ ∈ 0:ℓmax
    Vℓ′::Vector{Matrix{T}} # ℓ′ ∈ 1:ℓmax′
    reference::String
end
AtomicPotentials.charge(pp::PseudoPotential) = num_electrons(pp.gst_config)
AtomicPotentials.ground_state(pp::PseudoPotential) = pp.gst_config

Base.show(io::IO, pp::PseudoPotential) =
    write(io, "Relativistic pseudo-potential for $(pp.name) ($(pp.gst_config)), Z = $(charge(pp))")

function Base.show(io::IO, ::MIME"text/plain", pp::PseudoPotential)
    show(io, pp)
    println(io)
    ℓmax = length(pp.Vℓ)-1
    ℓmax′ = length(pp.Vℓ′)
    print(io, "Q = $(pp.Q), ℓ ∈ 0:$(ℓmax)")
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

function (pp::PseudoPotential{T})(orb::RelativisticOrbital, r::AbstractVector{T}) where T
    V = -pp.Q./r
    κ = orb.κ
    data = κ < 0 ? pp.Vℓ[-κ] : pp.Vℓ′[κ] # What to do if κ is too large?
    for k in 1:size(data,1)
        (n,β,B) = data[k,:]
        V += B * r.^(n-2).*exp.(-β*r.^2)
    end
    V
end
(pp::PseudoPotential{T})(orb::RelativisticOrbital, r::T) where T = pp([r])[1]

function parse_pseudopotential_line(::Type{T}, pp_line::AbstractString) where T
    data = split(pp_line, ";")
    num_terms = parse(Int, data[1])
    length(data) == num_terms+2 || throw(ArgumentError("Unrecognizable potential string $(pp_line)"))
    map(1:num_terms) do i
        term = split(strip(data[i+1]),",")
        parse.(Ref(T), term)'
    end |> p -> vcat(p...)
end

function parse_pseudopotential(::Type{T}, pp_string::String) where T
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

    Q = parse(Int, header[3])
    ncore = num_electrons(core(gst_config))
    Q == ncore ||
        throw(ArgumentError("Charge modelled by pseudopotential, Q = $Q, does not match number of core electrons in ground state configuration: $(core(gst_config)) => $(ncore) electrons"))

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
    PseudoPotential(element, gst_config, Q, Vℓ, Vℓ′, reference)
end

macro PP_str(pp_string)
    parse_pseudopotential(Float64, pp_string)
end

Neon = PP"""! [He] 2s2 2p6
!  Q=8., MEFIT, WB, Ref 22.
ECP,Ne,2,3,0;
1; 2,1.000000,0.00000000;
2; 2,31.860162,112.52543566; 2,12.362219,28.30083454;
2; 2,21.508034,-11.12658543; 2,12.910447,3.38754919;
1; 2,0.850385,-0.18408921;
! References:
! [22] A. Nicklass, M. Dolg, H. Stoll, H. Preuss, J. Chem. Phys. 102, 8942 (1995)."""

Argon = PP"""! [Ne] 3s2 3p6
!  Q=8., MEFIT, WB, Ref 22.
ECP,Ar,10,4,0;
1; 2,1.000000,0.00000000;
2; 2,10.261721,68.66778801; 2,3.952725,24.04276629;
2; 2,5.392714,27.73076331; 2,2.699967,4.04545904;
2; 2,8.086235,-8.13747696; 2,4.016632,-1.66452808;
1; 2,5.208459,-3.40009845;
! References:
! [22] A. Nicklass, M. Dolg, H. Stoll, H. Preuss, J. Chem. Phys. 102, 8942 (1995)."""

Krypton = PP"""! [Ne] 3s2 3p6 3d10 4s2 4p6
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

Xenon = PP"""! [Ar] 3d10c 4s2 4p6 4d10 5s2 5p6
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

Radon = PP"""! [Kr] 4d10c 4f14c 5s2 5p6 5d10 6s2 6p6
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

export @PP_str, @pc_str, charge, ground_state

end # module
