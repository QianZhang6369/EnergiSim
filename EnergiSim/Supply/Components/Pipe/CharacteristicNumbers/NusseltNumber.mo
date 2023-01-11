within EnergiSim.Supply.Components.Pipe.CharacteristicNumbers;
function NusseltNumber "Return Nusselt number"
  extends Modelica.Icons.Function;

  input Modelica.Units.SI.CoefficientOfHeatTransfer alpha "Coefficient of heat transfer";
  input Modelica.Units.SI.Length D "Characteristic dimension";
  input Modelica.Units.SI.ThermalConductivity lambda "Thermal conductivity";
  output Modelica.Units.SI.NusseltNumber Nu "Nusselt number";
algorithm
  Nu := alpha*D/lambda;
  annotation (Documentation(info="Nusselt number Nu = alpha*D/lambda"));
end NusseltNumber;
