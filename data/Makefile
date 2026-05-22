SATS=gps
test:
	cd .. && cargo build --release
	../target/release/tart-gnss-acquire --file data_2026-03-09_04_25_10.145497.hdf --${SATS} --filter-phase-mad 0.05
