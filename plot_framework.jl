using Plots

function plot_attributes!(adjust_legend=true)
    """
    This function gives some default attributes to plots. Run it after
        plotting anything you want to plot. It will customize the plot
        size to a golden ratio and adapt the fontsizes.
    """

    golden_ratio = (1+sqrt(5))/2
    width = 6.202*200               # Im Dokument mit so einem speziellen Paket gemessen. (Zoll)
    height = width/golden_ratio     # Sieht nice aus.

    plot!(dpi = 200, size = (width,height), guidefontsize = 25, tickfontsize = 15, legendfontsize = 20, titlefontsize = 30, widen = false,margin=15Plots.mm)
    adjust_legend ? plot!(legend = :topleft) : nothing
end
