within EnergiSim.Supply.Test;
model Pipe
  Components.Pipe.DynamicPipe pipe(
    redeclare package Medium = Modelica.Media.Water.ConstantPropertyLiquidWater,
    length=10,
    diameter=0.0127,
    nNodes=5,
    use_HeatTransfer=true,
    redeclare model HeatTransfer =
        EnergiSim.Supply.Components.Pipe.HeatTransfer.DBHeatTransfer)
    annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));

  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T=
        373.15)
    annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalCollector thermalCollector(m=
        5) annotation (Placement(transformation(extent={{-10,18},{10,-2}})));
  Modelica.Fluid.Sources.MassFlowSource_T boundary(
    redeclare package Medium = Modelica.Media.Water.ConstantPropertyLiquidWater,
    m_flow=0.1,
    T=303.15,
    nPorts=1)
    annotation (Placement(transformation(extent={{-66,-30},{-46,-10}})));

  Modelica.Fluid.Sources.FixedBoundary boundary1(redeclare package Medium =
        Modelica.Media.Water.ConstantPropertyLiquidWater, nPorts=1)
    annotation (Placement(transformation(extent={{64,-30},{44,-10}})));
equation
  connect(fixedTemperature.port, thermalCollector.port_b)
    annotation (Line(points={{-40,30},{0,30},{0,18}}, color={191,0,0}));
  connect(thermalCollector.port_a, pipe.heatPorts) annotation (Line(points={{0,
          -2},{0,-10.8},{0.1,-10.8},{0.1,-15.6}}, color={191,0,0}));
  connect(boundary.ports[1], pipe.port_a)
    annotation (Line(points={{-46,-20},{-10,-20}}, color={0,127,255}));
  connect(pipe.port_b, boundary1.ports[1])
    annotation (Line(points={{10,-20},{44,-20}}, color={0,127,255}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(StopTime=600, __Dymola_Algorithm="Dassl"));
end Pipe;
