within EnergiSim.Supply.Components.Pipe;
model DynamicPipe
  "Dynamic pipe model with storage of mass and energy"

  import Modelica.Fluid.Types.ModelStructure;

  // extending PartialStraightPipe
  extends EnergiSim.Supply.Components.Pipe.Basecase.PartialStraightPipe(
    final port_a_exposesState = (modelStructure == ModelStructure.av_b) or (modelStructure == ModelStructure.av_vb),
    final port_b_exposesState = (modelStructure == ModelStructure.a_vb) or (modelStructure == ModelStructure.av_vb));

  // extending PartialTwoPortFlow
  extends Modelica.Fluid.Pipes.BaseClasses.PartialTwoPortFlow(
    final lengths=fill(length/n, n),
    final crossAreas=fill(crossArea, n),
    final dimensions=fill(4*crossArea/perimeter, n),
    final roughnesses=fill(roughness, n),
    final dheights=height_ab*dxs);

  // Wall heat transfer
  parameter Boolean use_HeatTransfer = false
    "= true to use the HeatTransfer model"
      annotation (Dialog(tab="Assumptions", group="Heat transfer"));
  replaceable model HeatTransfer =
      Pipe.HeatTransfer.IdealFlowHeatTransfer
    constrainedby
    Pipe.HeatTransfer.PartialFlowHeatTransfer
    "Wall heat transfer"
      annotation (Dialog(tab="Assumptions", group="Heat transfer",enable=use_HeatTransfer),choicesAllMatching=true);
  Modelica.Fluid.Interfaces.HeatPorts_a[nNodes] heatPorts if use_HeatTransfer
    annotation (Placement(transformation(extent={{-10,45},{10,65}}),
        iconTransformation(extent={{-30,36},{32,52}})));

  HeatTransfer heatTransfer(
    redeclare final package Medium = Medium,
    final n=n,
    final nParallel=nParallel,
    final surfaceAreas=perimeter*lengths,
    final lengths=lengths,
    final dimensions=dimensions,
    final roughnesses=roughnesses,
    final states=mediums.state,
    final vs = vs,
    final use_k = use_HeatTransfer) "Heat transfer model"
      annotation (Placement(transformation(extent={{-45,20},{-23,42}})));
  final parameter Real[n] dxs = lengths/sum(lengths) "Normalized lengths";
equation
  Qb_flows = heatTransfer.Q_flows;
  // Wb_flow = v*A*dpdx + v*F_fric
  //         = v*A*dpdx + v*A*flowModel.dp_fg - v*A*dp_grav
  if n == 1 or useLumpedPressure then
    Wb_flows = dxs * ((vs*dxs)*(crossAreas*dxs)*((port_b.p - port_a.p) + sum(flowModel.dps_fg) - system.g*(dheights*mediums.d)))*nParallel;
  else
    if modelStructure == ModelStructure.av_vb or modelStructure == ModelStructure.av_b then
      Wb_flows[2:n-1] = {vs[i]*crossAreas[i]*((mediums[i+1].p - mediums[i-1].p)/2 + (flowModel.dps_fg[i-1]+flowModel.dps_fg[i])/2 - system.g*dheights[i]*mediums[i].d) for i in 2:n-1}*nParallel;
    else
      Wb_flows[2:n-1] = {vs[i]*crossAreas[i]*((mediums[i+1].p - mediums[i-1].p)/2 + (flowModel.dps_fg[i]+flowModel.dps_fg[i+1])/2 - system.g*dheights[i]*mediums[i].d) for i in 2:n-1}*nParallel;
    end if;
    if modelStructure == ModelStructure.av_vb then
      Wb_flows[1] = vs[1]*crossAreas[1]*((mediums[2].p - mediums[1].p)/2 + flowModel.dps_fg[1]/2 - system.g*dheights[1]*mediums[1].d)*nParallel;
      Wb_flows[n] = vs[n]*crossAreas[n]*((mediums[n].p - mediums[n-1].p)/2 + flowModel.dps_fg[n-1]/2 - system.g*dheights[n]*mediums[n].d)*nParallel;
    elseif modelStructure == ModelStructure.av_b then
      Wb_flows[1] = vs[1]*crossAreas[1]*((mediums[2].p - mediums[1].p)/2 + flowModel.dps_fg[1]/2 - system.g*dheights[1]*mediums[1].d)*nParallel;
      Wb_flows[n] = vs[n]*crossAreas[n]*((port_b.p - mediums[n-1].p)/1.5 + flowModel.dps_fg[n-1]/2+flowModel.dps_fg[n] - system.g*dheights[n]*mediums[n].d)*nParallel;
    elseif modelStructure == ModelStructure.a_vb then
      Wb_flows[1] = vs[1]*crossAreas[1]*((mediums[2].p - port_a.p)/1.5 + flowModel.dps_fg[1]+flowModel.dps_fg[2]/2 - system.g*dheights[1]*mediums[1].d)*nParallel;
      Wb_flows[n] = vs[n]*crossAreas[n]*((mediums[n].p - mediums[n-1].p)/2 + flowModel.dps_fg[n]/2 - system.g*dheights[n]*mediums[n].d)*nParallel;
    elseif modelStructure == ModelStructure.a_v_b then
      Wb_flows[1] = vs[1]*crossAreas[1]*((mediums[2].p - port_a.p)/1.5 + flowModel.dps_fg[1]+flowModel.dps_fg[2]/2 - system.g*dheights[1]*mediums[1].d)*nParallel;
      Wb_flows[n] = vs[n]*crossAreas[n]*((port_b.p - mediums[n-1].p)/1.5 + flowModel.dps_fg[n]/2+flowModel.dps_fg[n+1] - system.g*dheights[n]*mediums[n].d)*nParallel;
    else
      assert(false, "Unknown model structure");
    end if;
  end if;

  connect(heatPorts, heatTransfer.heatPorts)
    annotation (Line(points={{0,55},{0,54},{-34,54},{-34,38.7}}, color={191,0,0}));
  annotation (defaultComponentName="pipe",
Documentation(info="<html>
<p>Based on the model from Modelica.Fluid, implement the Dittus-Boelter heat transfer correlation</p>
<p>Model of a straight pipe with distributed mass, energy and momentum balances. It provides the complete balance equations for one-dimensional fluid flow as formulated in <a href=\"modelica://Modelica.Fluid.UsersGuide.ComponentDefinition.BalanceEquations\">UsersGuide.ComponentDefinition.BalanceEquations</a>.</p>
<p>This generic model offers a large number of combinations of possible parameter settings. In order to reduce model complexity, consider defining and/or using a tailored model for the application at hand, such as <a href=\"modelica://Modelica.Fluid.Examples.HeatExchanger.HeatExchangerSimulation\">HeatExchanger</a>.</p>
<p>DynamicPipe treats the partial differential equations with the finite volume method and a staggered grid scheme for momentum balances. The pipe is split into nNodes equally spaced segments along the flow path. The default value is nNodes=2. This results in two lumped mass and energy balances and one lumped momentum balance across the dynamic pipe.</p>
<p>Note that this generally leads to high-index DAEs for pressure states if dynamic pipes are directly connected to each other, or generally to models with storage exposing a thermodynamic state through the port. This may not be valid if the dynamic pipe is connected to a model with non-differentiable pressure, like a Sources.Boundary_pT with prescribed jumping pressure. The <b><span style=\"font-family: Courier New;\">modelStructure</span></b> can be configured as appropriate in such situations, in order to place a momentum balance between a pressure state of the pipe and a non-differentiable boundary condition.</p>
<p>The default <b><span style=\"font-family: Courier New;\">modelStructure</span></b> is <span style=\"font-family: Courier New;\">av_vb</span> (see Advanced tab). The simplest possible alternative symmetric configuration, avoiding potential high-index DAEs at the cost of the potential introduction of nonlinear equation systems, is obtained with the setting <span style=\"font-family: Courier New;\">nNodes=1, modelStructure=a_v_b</span>. Depending on the configured model structure, the first and the last pipe segment, or the flow path length of the first and the last momentum balance, are of half size. See the documentation of the base class <a href=\"modelica://Modelica.Fluid.Pipes.BaseClasses.PartialTwoPortFlow\">Pipes.BaseClasses.PartialTwoPortFlow</a>, also covering asymmetric configurations.</p>
<p>The <b><span style=\"font-family: Courier New;\">HeatTransfer</span></b> component specifies the source term <span style=\"font-family: Courier New;\">Qb_flows</span> of the energy balance. The default component uses a constant coefficient for the heat transfer between the bulk flow and the segment boundaries exposed through the <span style=\"font-family: Courier New;\">heatPorts</span>. The <span style=\"font-family: Courier New;\">HeatTransfer</span> model is replaceable and can be exchanged with any model extended from <a href=\"modelica://Modelica.Fluid.Pipes.BaseClasses.HeatTransfer.PartialFlowHeatTransfer\">BaseClasses.HeatTransfer.PartialFlowHeatTransfer</a>.</p>
<p>The intended use is for complex networks of pipes and other flow devices, like valves. See, e.g.,</p>
<ul>
<li><a href=\"modelica://Modelica.Fluid.Examples.BranchingDynamicPipes\">Examples.BranchingDynamicPipes</a>, or</li>
<li><a href=\"modelica://Modelica.Fluid.Examples.IncompressibleFluidNetwork\">Examples.IncompressibleFluidNetwork</a>.</li>
</ul>
</html>", revisions="<html>
<p>- 2023.01.11,by QianZhang_UiS</p>
<p>First implementation.</p>
</html>"),
Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{100,100}}), graphics={
        Rectangle(
          extent={{-100,44},{100,-44}},
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={0,127,255}),
        Ellipse(
          extent={{-72,10},{-52,-10}},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{50,10},{70,-10}},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-48,15},{46,-20}},
          textString="%nNodes")}),
Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{100,
            100}}), graphics={
        Rectangle(
          extent={{-100,60},{100,50}},
          fillColor={255,255,255},
          fillPattern=FillPattern.Backward),
        Rectangle(
          extent={{-100,-50},{100,-60}},
          fillColor={255,255,255},
          fillPattern=FillPattern.Backward),
        Line(
          points={{100,45},{100,50}},
          arrow={Arrow.None,Arrow.Filled},
          pattern=LinePattern.Dot),
        Line(
          points={{0,45},{0,50}},
          arrow={Arrow.None,Arrow.Filled},
          pattern=LinePattern.Dot),
        Line(
          points={{100,-45},{100,-50}},
          arrow={Arrow.None,Arrow.Filled},
          pattern=LinePattern.Dot),
        Line(
          points={{0,-45},{0,-50}},
          arrow={Arrow.None,Arrow.Filled},
          pattern=LinePattern.Dot),
        Line(
          points={{-50,60},{-50,50}},
          pattern=LinePattern.Dot),
        Line(
          points={{50,60},{50,50}},
          pattern=LinePattern.Dot),
        Line(
          points={{0,-50},{0,-60}},
          pattern=LinePattern.Dot)}));
end DynamicPipe;
