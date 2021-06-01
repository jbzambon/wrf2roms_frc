% D_WRF2ROMS:  Modified driver script from D_ECMWF2ROMS to create ROMS
%              forcing NetCDF from WRF output (via NCL script wrfout_to_cf.ncl)
% Joseph B. Zambon
%  jbzambon@ncsu.edu
%  14 May 2021
%
% D_ECMWF2ROMS:  Driver script to create ROMS forcing NetCDF from ECMWF
%                ERA-Interim Dataset
%
% This is a user modifiable script showing how to create ROMS forcing
% NetCDF files. The data source is the ECMWF's ERA-Interim Dataset
%
% http://apps.ecmwf.int/datasets/data/interim_full_daily/
%
% This global data is available from Jan 1, 1979 to the present. It is
% usually extracted to a regional grid from the ECMWF data server.  No
% attempt is made to interpolate the forcing fields to a particular
% ROMS grid.  ROMS has the capability to perform the spatial and
% temporal interpolation of 2D forcing fields.
%
% This is a template script showing how to create several NetCDF files
% (one per each atmospheric forcing field) using ROMS metadata structure,
% which follows the schema of "nc_inq" or native 'ncinfo" function.
%
% The following parameters are used to extract ERA-Interim fields:
%
% Select time:   00:00:00     12:00:00
%
% Select step:   0  3  6  9  12
%

% svn $Id: d_ecmwf2roms.m 796 2016-05-11 01:45:21Z arango $
%=========================================================================%
%  Copyright (c) 2002-2012 The ROMS/TOMS Group      Hernan G. Arango      %
%    Licensed under a MIT/X style license           John Wilkin           %
%    See License_ROMS.txt                                                 %
%=========================================================================%
%
% Set file names. Input ERA-Interim input data file names.  The data
% was extracted for the Gulf of Mexico (100.5W, 10.5N) to (70.5W, 40.5N) 
% every 3-hours. The dataset was extracted in several annual files:
%
%                                               ERA Variables
%
%   ecmwf_era_atms_2000.nc              time, msl, tcc, v10u, v10u
%   ecmwf_era_flux_2000.nc              time, ewss, nsss, e, tp
%   ecmwf_era_heat_2000.nc              time, sshf, slhf, ssr, str
%   ecmwf_era_temp_2000.nc              time, v2t, v2d, ssrd, par
%   ...
%   ecmwf_era_atms_2012.nc
%   ecmwf_era_flux_2012.nc
%   ecmwf_era_heat_2012.nc
%   ecmwf_era_temp_2012.nc
%
% ERA-Interim variable names:  @ denotes accumulated quantity that must
%                              be divided by the time interval into the
%                              cycle 3, 6, 9 or 12 hours
%
% @  sshf    W m-2 s        surface sensible heat flux
% @  slhf    W m-2 s        surface latent heat flux
% @  ssr     W m-2 s        surface net solar radiation (shortwave)
% @  str     W m-2 s        surface net thermal radiation (longwave)
% @  strd    W m-2 s        surface thermal radiation downwards
% @  ewss    N m-2 s        east-west surface stress
% @  nsss    N m-2 s        north-south surface stress
% @  e       m              evaporation (downward flux is positive)
% @  ro      m              runoff
% @  tcc     nondimensional total cloud cover [0:1]
% @  tp      m              total precipitation
% @  par     W m-2 s        photosynthetically active radiation at surface
%    msl     Pa             mean sea level pressure
%    v10v    m s-1          10 metre U wind component
%    vl0u    m s-1          10 metre V wind component
%    v2t     K              2 metre temperature
%    v2d     K              2 metre dewpoint temperature
%
% This dataset is written in compact way (short numbers). We need to
% convert to floating-point data and scale to ROMS units:
%
%   Uwind       (m s-1)         v10u
%   Vwind       (m s-1)         v10v
%   sustr       (N m-2)         ewss / (3*3600);   3-hour step
%   svstr       (N m-2)         nsss / (3*3600);   3-hour step
%   shflux      (W m-2)         (ssr+str+sshf+slhf) / (3*3600)
%   swrad       (W m-2)         ssr  / (3*3600);   3-hour step
%   lwrad_down  (W m-2)         strd / (3*3600);   3-hour step
%   latent      (W m-2)         slhf / (3*3600);   3-hour step
%   sensible    (W m-2)         sshf / (3*3600):   3-hour step
%   rain        (kg m-2 s-1)    tp * Rho_w / (3*3600)
%   evaporation (kg m-2 s-1)    e  * Rho_w / (3*3600)
%   swflux      (cm day-1)      (-e - tp) * 100 / (3/24);  0.125 day step
%   cloud       (nondimesional) tcc
%   Pair        (mb)            msl / 100;   (1 mb = 100 Pa)
%   Tair        (Celsius)       t2m - 273.15;   (1 C = 273.15 K)
%   Qair        (percentage)    100 * (E/Es)
%
% where
%
%   Rho_w = 1000 kg m-3  (water density)
%
%   E  = 6.11 * 10.0 ^ (7.5 * d2m / (237.7 + d2m))    vapor pressure (mb)
%                                                     d2m in Celsius
%
%   Es = 6.11 * 10.0 ^ (7.5 * t2m / (237.7 + t2m))    saturation vapor
%                                                     pressure (mb)
%                                                     t2m in Celsius
%
%--------------------------------------------------------------------------

% Path to downloaded ERA data files.

% Dir = fullfile(getenv('HOME'),                                          ...
%                'ocean/repository/Projects/useast/Data/Forcing2');

Dir = fullfile('/scratch/jbzambon/florence/oct20_exp/WRF/ndown/WRF/wrf_frc/');

fint = 1  %forecast intervals in hours (e.g. 1 = hourly output)
qair_calc_from_td = true    %calculate qair from dewpoint or use RH from NCL

% Base date for ROMS forcing files: "days since 1858-11-17 00:00:00".

mybasedate = datenum(1858,11,17,0,0,0);

% Set ROMS and ECMWF-ERA field NetCDF variable names for each forcing
% variable, input file name prefix, output file name. Notice that output
% forcing variables follows ROMS metadata design.
%
% If the ROMS option BULK_FLUXES is NOT activatived, we just need to
% process elements 1:5 in the structure vector (F).  Otherwise, we need
% process almost all the fields below.

clear S F

addpath('/scratch/jbzambon/florence/oct20_exp/roms_only/Tools/mfiles/rutgers/grid/')
addpath('/scratch/jbzambon/florence/oct20_exp/roms_only/Tools/mfiles/rutgers/utility/')
addpath('/scratch/jbzambon/florence/oct20_exp/roms_only/Tools/mfiles/rutgers/boundary/')
addpath('/scratch/jbzambon/florence/oct20_exp/roms_only/Projects/florence/downscale/mexcdf/mexnc')
addpath('/scratch/jbzambon/florence/oct20_exp/roms_only/Projects/florence/downscale/mexcdf/snctools')
addpath('/scratch/jbzambon/florence/oct20_exp/roms_only/Projects/florence/downscale/mexcdf/netcdf_toolbox/netcdf')
addpath('/scratch/jbzambon/florence/oct20_exp/roms_only/Projects/florence/downscale/mexcdf/netcdf_toolbox/netcdf/nctype')
addpath('/scratch/jbzambon/florence/oct20_exp/roms_only/Projects/florence/downscale/mexcdf/netcdf_toolbox/netcdf/ncutility')
addpath('/scratch/jbzambon/florence/oct20_exp/roms_only/Tools/mfiles/rutgers/netcdf/')
addpath('/scratch/jbzambon/florence/oct20_exp/roms_only/Tools/mfiles/rutgers/initial')

F(1).Vname  = {'Uwind', 'u_10m_tr'};
F(1).Tname  = {'wind_time', 'time'};
F(1).input  = 'init_post.nc';
F(1).output = 'useast_wind_era.nc';
F(1).scale  = 1.0;

F(2).Vname  = {'Vwind', 'v_10m_tr'};
F(2).Tname  = {'wind_time', 'time'};
F(2).input  = 'init_post.nc';
F(2).output = 'useast_wind_era.nc';
F(2).scale  = 1.0;

F(3).Vname  = {'swrad', 'ssr'};
F(3).Tname  = {'srf_time',  'time'};
F(3).input  = 'useast_ecmwf_';
F(3).output = 'useast_swrad_era.nc';
F(3).scale  = -1.0/(fint*3600.0);

F(4).Vname  = {'lwrad_down', 'strd'};
F(4).Tname  = {'lrf_time',  'time'};
F(4).input  = 'useast_ecmwf_';
F(4).output = 'useast_lwrad_era.nc';
F(4).scale  = -1.0/(fint*3600.0);

F(5).Vname  = {'rain', 'tp'};
F(5).Tname  = {'rain_time',  'time'};
F(5).input  = 'useast_ecmwf_';
F(5).output = 'useast_rain_era.nc';
F(5).scale  = -1000.0/(fint*3600.0);

F(6).Vname  = {'Pair', 'msl'};
F(6).Tname  = {'pair_time',  'time'};
F(6).input  = 'useast_ecmwf_';
F(6).output = 'useast_Pair_era.nc';
F(6).scale  = 0.01;

F(7).Vname  = {'Tair', 't2m'};
F(7).Tname  = {'tair_time',  'time'};
F(7).input  = 'useast_ecmwf_';
F(7).output = 'useast_Tair_era.nc';
F(7).scale  = 1.0;

F(8).Vname  = {'Qair', 'd2m'};         % Use temperature (v2t) and
F(8).Tname  = {'qair_time',  'time'};  % dewpoint temperature (v2d)
F(8).input  = 'useast_ecmwf_';         % to compute relative humidity
F(8).output = 'useast_Qair_era.nc';
F(8).scale  = 1.0;

% Set field element in structure F to process.

doFields = [1 2 3 4 5 6 7 8];

% Various parameters.

spherical = true;                       % Spherical switch

FlipLat   = false;                       % Flip data in ERA NetCDF files
                                        % (see below information)

StrDay    = datenum('11-Sept-2018');     % starting day to process
EndDay    = datenum('11-Sept-2018');     % ending   day to process

nctype    = 'nc_float';                 % Input data is in single precision
Unlimited = true;                       % time dimension is umlimited in
                                        % output files

%--------------------------------------------------------------------------
% Create surface forcing NetCDF files: build creation parameters
% structure, S.  We want to create a single file for each forcing
% field.  However, the wind components "Uwind" and "Vwind" are
% saved in the same NetCDF file.
%--------------------------------------------------------------------------

Tindex       = [];
ReplaceValue = NaN;
PreserveType = true;

mode = netcdf.getConstant('CLOBBER');                    % overwrite!!!
mode = bitor(mode,netcdf.getConstant('64BIT_OFFSET'));

% The strategy here is to build manually the NetCDF metadata structure to
% facilitate creating several NetCDF file in a generic and compact way.
% This structure is similar to that returned by "nc_inq" or native Matlab
% function "ncinfo".
%
% Notice that we call the "roms_metadata" function to create the fields
% in the structure.  Then, we call "check_metadata" for fill unassigned
% values and to check for consistency.

OutFiles = {};

for n = doFields,

  ncname = char(F(n).output);
  Vname  = char(F(n).Vname{1});
  Tname  = char(F(n).Tname{1});
  
  if (n == 1),
    create_file = true;
  else
    create_file = ~any(strcmp(OutFiles, ncname));
  end
  
  OutFiles = [OutFiles, ncname];
  InpFile  = fullfile(Dir,'init_post.nc'); % have only one file
  
  if (create_file),
    disp(blanks(1));
    disp(['** Creating ROMS NetCDF forcing file: ', ncname,' **']);
    disp(blanks(1))
  end
  
  S.Filename = ncname;

  lon = nc_read(InpFile, 'lon',                                   ...
                Tindex, ReplaceValue, PreserveType);

  lat = nc_read(InpFile, 'lat',                                    ...
                Tindex, ReplaceValue, PreserveType);

  ROMS_lon = lon;

  if (FlipLat),
    ROMS_lat = flipud(lat);
  else
    ROMS_lat = lat;
  end

  Im = length(ROMS_lon(1,:));
  Jm = length(ROMS_lat(:,1));

  S.Attributes(1).Name      = 'type';
  S.Attributes(1).Value     = 'FORCING file';

  S.Attributes(2).Name      = 'title';
  S.Attributes(2).Value     = ['NCSU USEast WRF'];

  S.Attributes(3).Name      = 'history';
  S.Attributes(3).Value     = ['Forcing file created with ',            ...
                               which(mfilename) ' on ' date_stamp];

  S.Dimensions(1).Name      = 'lon';
  S.Dimensions(1).Length    = Im;
  S.Dimensions(1).Unlimited = false;

  S.Dimensions(2).Name      = 'lat';
  S.Dimensions(2).Length    = Jm;
  S.Dimensions(2).Unlimited = false;

  S.Dimensions(3).Name      = Tname;
  S.Dimensions(3).Length    = nc_constant('nc_unlimited');
  S.Dimensions(3).Unlimited = true;

  S.Variables(1) = roms_metadata('spherical');
  S.Variables(2) = roms_metadata('lon');
  S.Variables(3) = roms_metadata('lat');
  S.Variables(4) = roms_metadata(Tname, [], [], Unlimited);
  S.Variables(5) = roms_metadata(Vname, spherical, nctype, Unlimited);

% Edit the time variable 'units" attribute for the correct reference
% time and add calendar attribute.

  natts = length(S.Variables(4).Attributes);
  iatt  = strcmp({S.Variables(4).Attributes.Name}, 'units');
  S.Variables(4).Attributes(iatt).Value = ['days since ' datestr(mybasedate,31)];

  S.Variables(4).Attributes(natts+1).Name  = 'calendar';
  S.Variables(4).Attributes(natts+1).Value = 'gregorian';

% Change the dimensions of "shflux", "swflux", "sustr" or "svstr" to
% (lon,lat) dimensions to allow interpolation from a coarser grid in
% ROMS. As default, 'roms_metadata' uses ROMS native coordinates.
%
% This illustrates also how to manipulate the default metadata values set
% in 'roms_metadata'.  It has to be done before calling 'check_metadata'.

  if (strcmp(Vname, 'shflux') ||                                        ...
      strcmp(Vname, 'swflux') ||                                        ...
      strcmp(Vname, 'sustr')  ||                                        ...
      strcmp(Vname, 'svstr')),
    S.Variables(5).Dimensions(1).Name = 'lon';
    S.Variables(5).Dimensions(2).Name = 'lat';
  end
  
% Check ROMS metadata structure.  Fill unassigned fields.

  S = check_metadata(S);
  
% Create forcing NetCDF files.  Write our coordinates. Notice that
% "nc_append" is used here since we wand both wind components "Uwind"
% and "Vwind" to be in the same output NetCDF file for CF conventions.

  if (create_file),
    ncid = nc_create(S.Filename, mode, S); % create a new NetCDF file

    %lon = repmat(ROMS_lon,  [1 Jm]);
    lon = ROMS_lon;
    %lat = repmat(ROMS_lat', [Im 1]);
    lat = ROMS_lat;
  
    status = nc_write(S.Filename, 'spherical', int32(spherical));
    status = nc_write(S.Filename, 'lon',       lon);
    status = nc_write(S.Filename, 'lat',       lat);
  else
    nc_append(S.Filename, S)               % append to existing NetCDF file
  end

end

%--------------------------------------------------------------------------
% Convert forcing data to floating-point and scale to ROMS units.
%--------------------------------------------------------------------------

% This particular data set has a time coordinate in days starting
% on 1-Jan-1900.

epoch = datenum('01-Jan-1900');

% Set reference time for this application using the specified base date.

ref_time = (mybasedate - epoch);

StrYear  = str2num(datestr(StrDay, 'yyyy'));
EndYear  = str2num(datestr(EndDay, 'yyyy'));

MyRec = zeros([1 length(F)]);

disp(blanks(1));

for n = doFields,

  year  = StrYear;

  while (year <= EndYear)
    InpFile = fullfile(Dir, strcat(char(F(n).input)));
    OutFile = char(F(n).output);
    Troms   = char(F(n).Tname{1});
    Tecmwf  = char(F(n).Tname{2});
    Vroms   = char(F(n).Vname{1});
    Vecmwf  = char(F(n).Vname{2});
    scale   = abs(F(n).scale);

% Determine time records to process.

    time = nc_read(InpFile, Tecmwf);
    time = time ./ 24;                              % hours to day

    StrRec = 1;
    EndRec = length(time);

% Read and write forcing fields. The forcing field are scaled to ROMS
% units.

    for Rec=StrRec:EndRec,

      MyRec(n) = MyRec(n) + 1;

      mydate = datestr(epoch+time(Rec));

      %disp(['** Processing: ', Vroms, '  for  ',mydate,' **']);  

      frc_time = time(Rec) - ref_time;

      switch Vroms
        case 'Tair'
          field = nc_read(InpFile, Vecmwf, Rec);
          field = field - 273.15;                   % Kelvin to Celsius
        case 'Qair'     
          if (qair_calc_from_td)
            tsur  = nc_read(InpFile, 'T_2m', Rec);     % 2m temperature
            tdew  = nc_read(InpFile, 'Td_2m', Rec);    % 2m dewpoint
            E     = 6.11 .* 10.0 .^ (7.5 .* tdew ./ (237.7 + tdew));
            Es    = 6.11 .* 10.0 .^ (7.5 .* tsur ./ (237.7 + tsur));
            field = 100.0 .* (E ./ Es);
          else
            field = nc_read(InpFile, 'rh_2m', Rec);
          end
        case 'swflux'     
          evap  = nc_read(InpFile, 'WaterFlx' , Rec);      % evaporation
          if rec==1
            prec = 0
          else
           rain_g = nc_read(InpFile, 'precip_g' , Rec) - nc_read(InpFile, 'precip_g' , Rec-1);
           rain_c = nc_read(InpFile, 'precip_c' , Rec) - nc_read(InpFile, 'precip_c' , Rec-1);
           prec   = (rain_g + rain_c) .* -1000.0/(fint*3600.0);      % precipitation
          end
          field = (-evap - prec) .* -1000.0/(fint*3600.0);
        case 'shflux'   
          sensbl = nc_read(InpFile, 'sshf' , Rec);  % sensible
          latent = nc_read(InpFile, 'slhf', Rec);   % latent
          nlwrad = nc_read(InpFile, 'LW_d', Rec);    % net longwave
          nswrad = nc_read(InpFile, 'SW_d', Rec);    % net shortwave
          field = (sensbl+latent+nlwrad+nswrad) .* scale;
        case 'rain'   %WRF does cumuluative rainfall, need to subtract previous timestep
          if rec==1
            rain_g = 0; rain_c=0;  field = (rain_g + rain_c) .* scale;   %rain at init =0
          else
            rain_g = nc_read(InpFile, 'precip_g' , Rec) - nc_read(InpFile, 'precip_g' , Rec-1);
            rain_c = nc_read(InpFile, 'precip_c' , Rec) - nc_read(InpFile, 'precip_c' , Rec-1);
            field = (rain_g + rain_c) .* scale;
          end
        otherwise
          field = nc_read(InpFile, Vecmwf, Rec);
          field = field.*scale;
      end

      fieldfinal = field;

% If the scale F(n).scale is set to negative, the input ECMWF data is a
% cummulative integral in forecast cycle from hour zero.
% For steps at 6, 9 and 12 hours we must separate last 3 hours of 
% integration from previous accumulation.
% At 3 hour step don't change anything

      if (F(n).scale < 0),
        step = rem(frc_time,0.5)*24;
        if step == 3
          fieldfinal = field;
        else
          fieldfinal = field - field_previous;  % At other steps subtract
        end                                     % the previous accumulation

        frc_time = frc_time - 1.5/24;	        % Center forcing time on the
                                                % accumulation interval

        field_previous = field;                 % Save this accumulation
                                                % to on the next step
      end

% If appropriate, flip the data so the origin (1,1) corresponds to
% the lower-left corner (LonMin,LatMin) of the extracted region.
% ERA data is written into NetCDF files with the origin at (1,end).
% That is, the origin is the left-upper corner of the extracted grid
% (LonMin,LatMax).  Users need to check and plot the extracted data
% from the ECMWF server.

      if (FlipLat),
        fieldfinal  = fliplr(fieldfinal);
      end

% Write out data.

      status = nc_write(OutFile, Troms, frc_time, MyRec(n));
      status = nc_write(OutFile, Vroms, fieldfinal, MyRec(n));
    end

%  Advance year counter.  Recall all input files are split in annual
%  files.

    year = year + 1;

  end
end
