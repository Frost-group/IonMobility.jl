using IonMobility
using Gnuplot

chi_values = collect(0.01:0.05:1.5)
D_ratios = [minimal_chakraborty_scpt(chi) for chi in chi_values]

@gp chi_values D_ratios "w l t 'D_{eff}/D'"
@gp :- "set xlabel 'λ_D * κ'"
@gp :- "set ylabel 'D_{eff}/D'"
@gp :- "set title 'Effective Diffusivity vs Disorder Strength'"
@gp :- "set grid"
@gp :- "set key top right"


Gnuplot.save("Chakraborty-Fig4-Deff.png", term="pngcairo size 800,600 enhanced font 'Helvetica,14'")
println("Plot saved to Chakraborty-Fig4-Deff.png")

Gnuplot.save("Chakraborty-Fig4-Deff.pdf", term="pdfcairo size 3in,2in enhanced font 'Helvetica,9'")
println("Plot saved to Chakraborty-Fig4-Deff.pdf")


