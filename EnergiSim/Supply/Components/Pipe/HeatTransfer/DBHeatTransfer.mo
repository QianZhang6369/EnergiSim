within EnergiSim.Supply.Components.Pipe.HeatTransfer;
model DBHeatTransfer
  extends PartialPipeFlowHeatTransfer;

equation

  for i in 1:n loop
   Nus[i]=0.023 * Prs[i] ^ (0.4) * Res[i] ^ (4 / 5);
  end for;

end DBHeatTransfer;
