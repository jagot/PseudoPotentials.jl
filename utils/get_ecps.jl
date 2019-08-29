using Pkg
Pkg.activate(".")

using AtomicLevels
using AtomicPotentials

using HTTP
using Gumbo
using AbstractTrees

const BASE_URL = "http://www.tc.uni-koeln.de/cgi-bin/pp.pl"

function ask_cologne(;kwargs...)
    kwargs = Dict(kwargs...)
    kwargs[:format] = get(kwargs, :format, "molpro")
    kwargs[:language] = get(kwargs, :language, "en")

    url = "$(BASE_URL)?"*join(["$(k)=$(w)" for (k,w) in kwargs], ",")
    r = HTTP.request("GET", url)
    if r.status == 200
        parsehtml(String(r.body))
    else
        throw(ArgumentError("Unexpected HTTP status: $(r.status)"))
    end
end

function get_ecps_list(element::String)
    doc = ask_cologne(job="listecps", element=element)

    ecps = String[]
    for elem in PreOrderDFS(doc.root)
        if elem isa HTMLElement{:a}
            a = attrs(elem)
            if "href" in keys(a) && occursin("getecp",a["href"])
                push!(ecps, first(Gumbo.children(elem)).text)
            end
        end
    end
    ecps
end

function core_model(element::String, Q::Integer)
    Z = element_number(Symbol(element))
    gst = ground_state(PointCharge(Z))
    o = 0
    for i in 1:length(gst)
        o += gst.occupancy[i]
        gst.states[i] = Z - o < Q ? :open : :closed
    end
    gst
end

function get_ecp(element::String, ecp::String)
    doc = ask_cologne(job="getecp", element=element, ecp=ecp)
    data = String[]
    for elem in PreOrderDFS(doc.root)
        if elem isa HTMLElement{:pre}
            push!(data, strip(first(Gumbo.children(elem)).text))
        end
    end
    m = match(r"Q=(.*?),", data[1])
    if !isnothing(m)
        Q = parse(Float64, m[1])
        if isinteger(Q)
            c = core_model(element, Int(Q))
            prepend!(data, ["! $(ascii(c))"])
        else
            @warn "Can't generate model configuration for non-integer Q = $Q"
        end
    end
    data = join(data, "\n")
    data
end

function get_ecps(element::String)
    ecps = get_ecps_list(element)
    if isempty(ecps)
        println("No ECPs available for $(element)")
        return
    end
    println(repeat("=", 100))
    println("""ECPnXY
n = number of core electrons
X = S/M: single/multi electron fit
Y = HF/WB/DF: non/quasi/fully relativistic""")
    for ecp in ecps
        println(repeat("=", 100))
        println("==> $(ecp)")
        println(get_ecp(element, ecp))
    end
    println(repeat("=", 100))
end
