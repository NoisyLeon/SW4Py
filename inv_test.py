import obspy

net=obspy.core.inventory.network.Network('SW4', stations=[])
inv=obspy.core.inventory.inventory.Inventory(networks=[net],source='CU')
sta=obspy.core.inventory.station.Station('aa',13,132.4214,0.0)
