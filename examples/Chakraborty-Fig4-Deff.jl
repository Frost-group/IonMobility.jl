using IonMobility
using Gnuplot

# function minimal_chakraborty_scpt(
#   λκ::Float64;
#   N::Int = 10,
#   T::Float64 = 10.0,
#   maxiter::Int = 200, )

@gp "set title 'Effective Diffusivity vs Disorder Strength'"
@gp :- "set xlabel 'λ_D * κ'"
@gp :- "set ylabel 'D_{eff}/D'"
@gp :- "set grid"
@gp :- "set key top right"

chi_values = collect(0.01:0.05:3.0)

for T in [10.0,25.0,50.0,100.0]
    D_ratios = [minimal_chakraborty_scpt(chi, T=T) for chi in chi_values]
    @gp :- chi_values D_ratios "w lp t 'T: $T D_{eff}/D'"
end

Gnuplot.save("Chakraborty-Fig4-Deff.png", term="pngcairo size 800,600 enhanced font 'Helvetica,14'")
println("Plot saved to Chakraborty-Fig4-Deff.png")

Gnuplot.save("Chakraborty-Fig4-Deff.pdf", term="pdfcairo size 3in,2in enhanced font 'Helvetica,9'")
println("Plot saved to Chakraborty-Fig4-Deff.pdf")


