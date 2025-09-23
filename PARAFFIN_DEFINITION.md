I had a bug I couldn't figure out how to fix with adding paraffin to rocketcea:
add_new_fuel(fuel_name, fuel_properties)

To overcome this issue, I added the following card to "input_cards.py" in rocketcea.

The information was obtained from NASA-SP-1311-2

#NOTE: RW MODIFICATION 21-SEPT-2025
fuelCards["paraffin"] = [" fuel paraffin C 1 H 2   wt%=100.00",
                         " h,cal=-25600.0   t(k)=298.15   rho=0.900 "]
