using DataFrames
using CSV
using Plots
using Dates
include("plot_framework.jl")

function transpose_df(dataframe)
    dataframe = DataFrame([[names(dataframe)]; collect.(eachrow(dataframe))], [:column; Symbol.(axes(dataframe, 1))])
end

# read in the hopkins data
dfconfirmed = DataFrame(CSV.file("../data_hopkins/time_series_19-covid-Confirmed.csv"))
dfdead = DataFrame(CSV.file("../data_hopkins/time_series_19-covid-Deaths.csv"))
dfrecovered = DataFrame(CSV.file("../data_hopkins/time_series_19-covid-Recovered.csv"))

# Storing country and region states
cnames = dfrecovered[Symbol("Country/Region")]
rnames = dfrecovered[Symbol("Province/State")]

# Every province/state value missing will be replaced by country
for (c,i) ∈ zip(rnames,1:length(rnames))
    if c === missing
        rnames[i] = cnames[i]
    end
end

# write the list of provinces/countrys into a columnname list
colnames = [:date]
for elem in rnames
    push!(colnames,Symbol(elem))
end

# transpose the dataframes for better handling
dfrecovered = transpose_df(dfrecovered)
dfconfirmed = transpose_df(dfconfirmed)
dfdead = transpose_df(dfdead)

# kick out unwanted information
deleterows!(dfconfirmed,1:4)
deleterows!(dfdead,1:4)
deleterows!(dfrecovered,1:4)

# renaming the rows
rename!(dfconfirmed, colnames,makeunique=true)
rename!(dfdead, colnames,makeunique=true)
rename!(dfrecovered, colnames, makeunique=true)

# plotting the data for each province
for col in filter!(e->e≠:date,names(dfconfirmed))
    plot()
    for (frame,ccd) in zip([dfconfirmed,dfdead,dfrecovered],["confirmed", "dead", "cured"])
        dfred = DataFrame(date = frame.date, cases = frame[:,col])
        dfred = dropmissing(dfred)
        a = split.(string.(dfred.date), "/")
        c = deepcopy(a)
        for (b,i) ∈ zip(c,1:length(c))
            a[i][1] = "20"*b[3]
            a[i][2] = b[1]
            a[i][3] = b[2]
        end

        dates = []

        for k in 1:length(a)
            b = a[k]
            push!(dates,b[1]*"-"*b[2]*"-"*b[3])
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
        CSV.write("../data_hopkins/$col-$ccd.csv",dfred)
        plot!(title = "# of cases in $col",yaxis="population")
        plot_attributes!()
        savefig("../data_hopkins/$col.pdf")
    end
end
