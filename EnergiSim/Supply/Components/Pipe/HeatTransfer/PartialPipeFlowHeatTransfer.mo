within EnergiSim.Supply.Components.Pipe.HeatTransfer;
partial model PartialPipeFlowHeatTransfer
    "Base class for pipe heat transfer correlation in terms of Nusselt number heat transfer in a circular pipe for laminar and turbulent one-phase flow"
  extends PartialFlowHeatTransfer;
  parameter Modelica.Units.SI.CoefficientOfHeatTransfer alpha0=100
      "Guess value for heat transfer coefficients";
  Modelica.Units.SI.CoefficientOfHeatTransfer[n] alphas(each start=alpha0)
      "Heat transfer coefficient";
  Real[n] Res "Reynolds numbers";
  Real[n] Prs "Prandtl numbers";
  Real[n] Nus "Nusselt numbers";
  Medium.Density[n] ds "Densities";
  Medium.DynamicViscosity[n] mus "Dynamic viscosities";
  Medium.ThermalConductivity[n] lambdas "Thermal conductivity";
  Modelica.Units.SI.Length[n] diameters = dimensions "Hydraulic diameters for pipe flow";
equation
  ds=Medium.density(states);
  mus=Medium.dynamicViscosity(states);
  lambdas=Medium.thermalConductivity(states);
  Prs = Medium.prandtlNumber(states);
  Res =
    Pipe.CharacteristicNumbers.ReynoldsNumber(
    vs,
    ds,
    mus,
    diameters);
  Nus =
    Pipe.CharacteristicNumbers.NusseltNumber(
    alphas,
    diameters,
    lambdas);
  Q_flows={alphas[i]*surfaceAreas[i]*(heatPorts[i].T - Ts[i])*nParallel for i in 1:n};
    annotation (Documentation(info="<html>
<p>
Base class for heat transfer models that are expressed in terms of the Nusselt number and which can be used in distributed pipe models.
</p>
</html>"));
end PartialPipeFlowHeatTransfer;
