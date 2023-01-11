within EnergiSim.Supply.Components.Pipe.CharacteristicNumbers;
function ReynoldsNumber_m_flow
  "Return Reynolds number from m_flow, mu, D, A"
  extends Modelica.Icons.Function;

  input Modelica.Units.SI.MassFlowRate m_flow "Mass flow rate";
  input Modelica.Units.SI.DynamicViscosity mu "Dynamic viscosity";
  input Modelica.Units.SI.Length D
    "Characteristic dimension (hydraulic diameter of pipes or orifices)";
  input Modelica.Units.SI.Area A = Modelica.Constants.pi/4*D*D
    "Cross sectional area of fluid flow";
  output Modelica.Units.SI.ReynoldsNumber Re "Reynolds number";
algorithm
  Re := abs(m_flow)*D/A/mu;
  annotation (Documentation(info="<html>Simplified calculation of Reynolds Number for flow through pipes or orifices;
              using the mass flow rate <code>m_flow</code> instead of the velocity <code>v</code> to express inertial forces.
<blockquote><pre>
  Re = |m_flow|*diameter/A/&mu;
with
  m_flow = v*&rho;*A
</pre></blockquote>
See also <a href=\"modelica://Modelica.Fluid.Pipes.BaseClasses.CharacteristicNumbers.ReynoldsNumber\">
          Pipes.BaseClasses.CharacteristicNumbers.ReynoldsNumber</a>.
</html>"));
end ReynoldsNumber_m_flow;
