% time integration for LDR

halflife = machine.data.SourceIsotopeHalfLife;
TimeIntegrationFactor = 24*halflife/log(2)*(1-exp(-t*log(2)/halflife));