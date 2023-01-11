within EnergiSim.Supply.Components.Pipe.CharacteristicNumbers;
function ReynoldsNumber "Return Reynolds number from v, rho, mu, D"
  extends Modelica.Icons.Function;

  input Modelica.Units.SI.Velocity v "Mean velocity of fluid flow";
  input Modelica.Units.SI.Density rho "Fluid density";
  input Modelica.Units.SI.DynamicViscosity mu "Dynamic (absolute) viscosity";
  input Modelica.Units.SI.Length D
    "Characteristic dimension (hydraulic diameter of pipes)";
  output Modelica.Units.SI.ReynoldsNumber Re "Reynolds number";
algorithm
  Re := abs(v)*rho*D/mu;
  annotation (Documentation(info="<html>
<p>
Calculation of Reynolds Number
</p>
<blockquote><pre>
Re = |v|&rho;D/&mu;
</pre></blockquote>
<p>
a measure of the relationship between inertial forces (v&rho;) and viscous forces (D/&mu;).
</p>
<p>
The following table gives examples for the characteristic dimension D and the velocity v for different fluid flow devices:
</p>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th><strong>Device Type</strong></th><th><strong>Characteristic Dimension D</strong></th><th><strong>Velocity v</strong></th></tr>
<tr><td>Circular Pipe</td><td>diameter</td>
    <td>m_flow/&rho;/crossArea</td></tr>
<tr><td>Rectangular Duct</td><td>4*crossArea/perimeter</td>
    <td>m_flow/&rho;/crossArea</td></tr>
<tr><td>Wide Duct</td><td>distance between narrow, parallel walls</td>
    <td>m_flow/&rho;/crossArea</td></tr>
<tr><td>Packed Bed</td><td>diameterOfSpericalParticles/(1-fluidFractionOfTotalVolume)</td>
    <td>m_flow/&rho;/crossArea (without particles)</td></tr>
<tr><td>Device with rotating agitator</td><td>diameterOfRotor</td>
    <td>RotationalSpeed*diameterOfRotor</td></tr>
</table>
</html>"));
end ReynoldsNumber;
