SATS=gps
test:
	cd .. && RUSTFLAGS="-C target-cpu=native" cargo build --release
	../target/release/tart-gnss-acquire --file data_2026-03-09_04_25_10.145497.hdf --${SATS} --ant 0,1,2,3,4 --filter-phase-mad 0.05
