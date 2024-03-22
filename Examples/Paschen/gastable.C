#include "Garfield/MediumMagboltz.hh"

using namespace Garfield;

int main(int argc, char *argv[]) {
  MediumMagboltz* gas = new MediumMagboltz();
  gas->SetComposition("ar", 90.0, "ic4h10", 10.0);
  gas->SetTemperature(293.15);
  gas->SetPressure(760.);
  gas->SetFieldGrid(1000., 1500000., 40, true, 0., 0., 1, 0., 0., 1);
  const double lP = 0.0;
  const double rP = 0.0;
  gas->EnablePenningTransfer(rP, lP, "ar");
  gas->GenerateGasTable(3);
	
  gas->WriteGasFile("ar-90-ic4h10-10-p1.gas");
  return 0;
}
