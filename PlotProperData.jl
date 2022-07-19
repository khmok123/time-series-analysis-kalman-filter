using DataFrames
using CSV
using Plots
using Dates
include("plot_framework.jl")


# load csv data
df = DataFrame(CSV.File("../data/DXYArea.csv"))

# make new empty dataframes
confirmed = DataFrame()
dead = DataFrame()
cured = DataFrame()

# get province names and unique update times
provinces = unique(df.provinceEnglishName)
updateTimes = unique(df.updateTime)

# make new lists for each day to put into new dataframes
for (i,t) in enumerate(updateTimes)
    lstconfirmed = []
    lstdead = []
    lstcured = []
    for prov in provinces
        fdf = df[(df.provinceEnglishName .== prov) .& (df.updateTime .== t), :]
        confirmedcount = fdf.province_confirmedCount
        deadcount = fdf.province_deadCount
        curedcount = fdf.province_curedCount
        if length(confirmedcount) > 0
            push!(lstconfirmed,confirmedcount[1])
        else
            push!(lstconfirmed,missing)
        end
        if length(deadcount) > 0
            push!(lstdead,deadcount[1])
        else
            push!(lstdead,missing)
        end
        if length(curedcount) > 0
            push!(lstcured,curedcount[1])
        else
            push!(lstcured,missing)
        end
    end
    confirmed.time = lstconfirmed
    cured.time = lstcured
    dead.time = lstdead
    rename!(confirmed, Symbol.(updateTimes[1:i]))
    rename!(dead, Symbol.(updateTimes[1:i]))
    rename!(cured, Symbol.(updateTimes[1:i]))
end

# Transpose dataframes
cured = DataFrame([[names(cured)]; collect.(eachrow(cured))], [:column; Symbol.(axes(cured, 1))])
confirmed = DataFrame([[names(confirmed)]; collect.(eachrow(confirmed))], [:column; Symbol.(axes(confirmed, 1))])
dead = DataFrame([[names(dead)]; collect.(eachrow(dead))], [:column; Symbol.(axes(dead, 1))])

# create new column names
colnames = [Symbol("date")]
for elem in provinces
    push!(colnames,Symbol(elem))
end

# change column names
rename!(confirmed, colnames)
rename!(dead, colnames)
rename!(cured, colnames)

# Plotting the data for each province
for col in filter!(e->e≠:date,names(confirmed))
    plot()
    for (frame,ccd) in zip([confirmed,dead,cured],["confirmed", "dead", "cured"])
        dfred = DataFrame(date = frame.date, cases = frame[:,col])
        dfred = dropmissing(dfred)
        a = split.(string.(dfred.date))

        dates = []

        for k in 1:length(a)
            b = a[k]
            push!(dates,b[1])
        end

        dates = DateTime.(dates)

        dfred.date = dates

        # delete rows with multiple entries for the same day, choose maximum of that day
        dfred = by(dfred, :date) do sbdf
            max = maximum(sbdf.cases)
            (size(sbdf, 1)>1) ? sbdf[sbdf.cases.==max,:] : sbdf
        end
        # delete rows with same entries for same day (duplicates)
        unique!(dfred)
        scatter!(dfred.date, dfred.cases, label = "$ccd", markersize = 5, markershape = :d)
        CSV.write("../data/$col-$ccd.csv",dfred)
        plot!(title = "# of cases in $col")
        plot_attributes!()
        savefig("../data/$col.pdf")
    end
end
