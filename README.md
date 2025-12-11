# mjo-tracking

A modular MATLAB toolkit for identifying, tracking, and compositing Madden–Julian Oscillation (MJO) events using daily OLR data.  
This repository provides a clean workflow from preprocessing to event detection, propagation-speed estimation, and composite analysis.

---

## Folder Structure

src/
│
├── preprocess_olr.m % Preprocess OLR: trim years, remove leap days, save new NetCDF
├── compute_equatorial_olr.m % Compute ±5° / ±10° equatorial mean OLR
├── build_mjo_segments.m % Build Hovmöller slices and identify t0 (MJO minima)
├── track_mjo_all_years.m % Core algorithm: fit propagation lines & compute MJO phase speeds
│
├── classify_mjo_speed.m % Classify events into slow / fast groups
├── composite_olr.m % Compute OLR composites (±30 days)
├── ttest_composite.m % Student-t significance testing for composites
│
├── plot_phase_speed_hist.m % Histogram of MJO phase speeds
└── plot_mjo_hovmoller.m % WK99-style longitude–time composite plotting


---

## Requirements

- MATLAB R2020a or later  
- Statistics & Machine Learning Toolbox  
- NetCDF support (built-in for modern MATLAB)

---

## Input Data Format

The raw dataset should be a daily OLR NetCDF file, typically from NOAA:

olr(lat, lon, time)
lat(lat)
lon(lon)
time(time) % hours or days since a reference epoch


Place your raw file under:

data/raw/


---

## Typical Workflow

Below is a minimal example showing how to run the full pipeline.

```matlab
%% --------------------------------------------------------------
% 1. Preprocess OLR (remove leap days, extract 1979–2013)
%% --------------------------------------------------------------
preprocess_olr('data/raw/olr.day.mean.nc', ...
               'data/processed/olr_1979_2013.nc', ...
               1979, 2013);

%% --------------------------------------------------------------
% 2. Compute equatorial OLR
%% --------------------------------------------------------------
[olr_EQ, time, lon] = compute_equatorial_olr('data/processed/olr_1979_2013.nc');

%% --------------------------------------------------------------
% 3. Build Hovmöller segments and identify t0 events
%% --------------------------------------------------------------
Seg = build_mjo_segments(olr_EQ, time, lon);

%% --------------------------------------------------------------
% 4. Track MJO propagation for all years
%% --------------------------------------------------------------
K_real = 1:0.1:25;
ref_lon = 90;
data_all = track_mjo_all_years(Seg, K_real, ref_lon);

%% --------------------------------------------------------------
% 5. Classify fast vs slow MJO events
%% --------------------------------------------------------------
[data_all, data_fast, data_slow] = classify_mjo_speed( ...
    'PropagationInfo.xlsx', ...
    'PropagationInfo_all.xlsx', ...
    'PropagationInfo_fast.xlsx', ...
    'PropagationInfo_slow.xlsx');

%% --------------------------------------------------------------
% 6. Composite analysis
%% --------------------------------------------------------------
window = 30;
[AllComp, FastComp, SlowComp] = composite_olr( ...
    time, olr_EQ, data_all, data_fast, data_slow, window);

%% --------------------------------------------------------------
% 7. Student-t significance testing
%% --------------------------------------------------------------
[t_fast, t_slow] = ttest_composite(FastComp, SlowComp);

%% --------------------------------------------------------------
% 8. Plot figures
%% --------------------------------------------------------------
plot_phase_speed_hist(data_all(:,4));  % histogram of speeds
plot_mjo_hovmoller(FastComp, SlowComp, t_fast, t_slow, lon);

---

## Output

The pipeline produces:

Estimated MJO propagation speeds (m/s)

Classification of fast vs slow events

Composite OLR fields (all/fast/slow)

Significance masks (t-test)

Figures including:

Phase speed histogram

WK99-style OLR variance composites

Hovmöller diagrams

Figures can be placed under:

figures/

## License

MIT License.
