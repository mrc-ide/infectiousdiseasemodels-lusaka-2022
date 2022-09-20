## Example 1: Vaccination at birth

This is a simple simulation of a vaccine that is introduced at or soon after birth. Think back to your open SIR model. The aim of vaccinating at birth is to effectively remove new-born individuals from the \(Susceptible\) population before they come in contact with \(Infected\) individuals. So the ODE equations will be

$$\frac{dS}{dt} = (1 - p) bN -\beta \frac{SI}{N} +\mu S$$

$$\frac{dI}{dt} = \beta \frac{SI}{N} -\gamma I -\mu I$$

$$\frac{dR}{dt} = pbN +\gamma I -\mu R$$

Thus note that \(p\) represents the proportion of new-borns being vaccinated. Code in these ODEs into the editor and initialise it with the parameters provided for birth, death, recovery and transmission rates. Run the model for 200 days. Note we are analysing an additional model outputs, \(Reff\).
