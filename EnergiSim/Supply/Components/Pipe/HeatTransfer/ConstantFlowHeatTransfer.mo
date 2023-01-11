within EnergiSim.Supply.Components.Pipe.HeatTransfer;
model ConstantFlowHeatTransfer
    "ConstantHeatTransfer: Constant heat transfer coefficient"
  extends PartialFlowHeatTransfer;
  parameter Modelica.Units.SI.CoefficientOfHeatTransfer alpha0 "Heat transfer coefficient";
equation
  Q_flows = {alpha0*surfaceAreas[i]*(heatPorts[i].T - Ts[i])*nParallel for i in 1:n};
  annotation(Documentation(info="<html>
<p>
Simple heat transfer correlation with constant heat transfer coefficient, used as default component in distributed pipe models.
</p>
</html>"));
end ConstantFlowHeatTransfer;
