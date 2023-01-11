within EnergiSim.Supply.Components.Pipe.HeatTransfer;
model IdealFlowHeatTransfer
    "IdealHeatTransfer: Ideal heat transfer without thermal resistance"
  extends PartialFlowHeatTransfer;
equation
  Ts = heatPorts.T;
  annotation(Documentation(info="<html>
Ideal heat transfer without thermal resistance.
</html>"));
end IdealFlowHeatTransfer;
