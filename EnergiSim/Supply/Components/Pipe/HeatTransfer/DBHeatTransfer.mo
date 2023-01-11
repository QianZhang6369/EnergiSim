within EnergiSim.Supply.Components.Pipe.HeatTransfer;
model DBHeatTransfer
  extends PartialPipeFlowHeatTransfer;

equation

  for i in 1:n loop
   Nus[i]=0.023 * Prs[i] ^ (0.4) * Res[i] ^ (4 / 5);
  end for;

  annotation (Documentation(revisions="<html>
<p>- 2023.01.11,by QianZhang</p>
<p>First implementation.</p>
</html>", info="<html>
<p>This is the model for heat transfer correlation (Dittus-Boelter)</p>
</html>"));
end DBHeatTransfer;
