# Log file - Telescope operations

**DO NOT MODIFY THE SCTRUCTURE OF THIS FILE, AS IT IS HANDELED IN update_daily_fluxes**

This file has been started by Andrea F. Battaglia, with the goal of keeping track of manual or technical telescope operations that may affect the observations. If there is no entry for a particular date means that the nominal solar measurements are done.
**THIS DOCUMENT HAS BEEN CREATED ON 03/12/2024, therefore there is no track on particular telescope operations before this date.**
For any questions, please contact andrea.francesco.battaglia@usi.irsol.ch


## Summary table

| **Date and time (UT)** | **Operator**  | **Comments**                       |
| ---------------------- | ------------- | ---------------------------------- |
| 2024/12/03 12:50       | Andrea        | Calibration of the pointing offset (noon) |
| 2024/12/03 14:39       | Andrea        | Calibration of the pointing offset (afternoon). The flux (in the evening) increased of about 50-70%! However, we noticed that in the offsets did not change from noon to afternoon |
| 2024/12/04 09:50       | Andrea        | Calibration of the pointing offset (morning) |
| 2024/12/09 09:20       | Andrea        | Calibration of the pointing offset (morning) |
| 2024/12/09 11:33       | Andrea        | Calibration of the pointing offset (noon) |
| 2024/12/09 14:35       | Andrea        | Calibration of the pointing offset (afternoon) |
| 2024/12/17 14:00       | Andrea        | Starting now, we upload the frq99995.cfg file to the automatic procedure. This shows the correct frequencies. Before, we uploaded frq00001.cfg, which is wrong! |
| 2024/12/18 09:20       | Maurizio      | We had to switch off the LN (RF survey) |
| 2024/12/18 10:15       | Andrea        | Switch callisto off, for manual pointing (RF survey) |
| 2024/12/19 12:15       | Marco         | Cleaning of the dish |
| 2024/12/23 --:--       | Andrea        | Upload the scripts for fixed Thot and Tcold positions. Some tests have been run throughout the day. |
| 2025/01/27 09:50       | Andrea & Maurizio  | 1. Switch off/on CALLISTO; 2. Switch Callisto OFF, disconnect N Connector, Switch Callisto ON; Switch OFF Callisto, reconnect N Connector, Switch ON Callisto. From preliminary analysis, it seems that after this operation, the overall signal is increased by 1 digital unit, hence about 10 mV. This may imply that, after this operation, the antenna gain (compared to the measured performed eariler) may be different. This needs to be checked after about a month of observations. |
| 2025/02/12 08:00       | CALLISTO | CALLISTO crashed. It is back working on 2025/02/14. |
| 2025/02/25 08:00       | Andrea & Philipp | Estimation of the elevation gradient, before the start time of the monitoring measurements. |
| 2025/03/03 14:45       | Andrea | Calibration of the pointing offset |
| 2025/03/05 10:30       | Andrea & Marco | Interference with the beam path (assessment of the alignment) |
| 2025/03/05 14:30       | Andrea | Calibration of the pointing offset |
| 2025/03/10 10:45       | CALLISTO SW | The software crashed. The FIT files generated on this day are corrupted |
| 2025/03/11 07:30       | Andsrea | Completely rebooted PC and CALLISTO spectrometer. Everything seems to be back running |
| 2025/03/11 08:00       | Andrea | Unfolded and folded the motor protection, for taking pictures as requested by Philipp. This day, lower fluxes have been recorded: could it be due to this? |
| 2025/03/11 09:45       | Andrea | Calibration of the pointing offset |
| 2025/03/11 15:00       | Andrea | After the calibration of the pointing offset at 09:45, the data look really odd. Because of this, I've reset the offset values in place before 09:45. |
| 2025/03/14 14:00       | Andrea | <span style="color:red"> **IMPORTANT:** solar data from 2025/03/11 to 2025/03/14 (included) are useless, because after the calibration of the pointing offset I have not set back CaliProc=True. It has been stayed CaliProc=False for three consecutive days! </span> |
| 2025/03/14 14:15       | Andrea | Because of the issue stated just here above, I have back the pointing offset as it was measured on 2025/03/11 at 09:45 |
| 2025/03/17 09:34       | Andrea | In the configsun file I set CaliProc = False and scanning = True. This is done in order to find out the best offset value for the calibration of the pointing. |
| 2025/03/17 14:50       | Andrea | Set back CaliProc = True and scanning = False. |
| 2025/03/18 08:45       | Andrea | Due to the scanning test of yesterday, I adapted aziref and eleref. I want to check now if the solar flux improved. |
| 2025/03/18 09:45       | Andrea | Changed the aziref further. |
| 2025/03/18 10:04       | Andrea | In the configsun file I set CaliProc = False and scanning = True. This is done in order to find out the best offset value for the calibration of the pointing. Test-day number 2. |
| 2025/03/18 14:45       | Andrea | Implementation of the scanning59 mode in sunpos_AZI_ELE.py. At 08:59, 10:59 and 12:59 CALLISTO does a scan of the sun to assess the pointing of the telescope. |
| 2025/03/19 08:55       | Andrea | Manual calibration of the pointing offset: aziref and eleref have been changed. |
| 2025/03/20 07:23       | Andrea | Slight changes to the elevations of Tcold and Thot, now to 72 deg (above summer solstice) and -15 deg (further down), respectively. |
| 2025/03/20 12:30       | Andrea | The elevation at -15 deg for Thot was out of range! In the morning, no Thot measurements have been taken. Now it has been changed to -14 deg (the minimum elevation is -15 deg). |
| 2025/03/21 08:00       | Andrea, Maurizio & Phil | Test of the CALLISTO spectrometer. We run measurements until 10:20 UT. |
| 2025/03/21 13:50       | Andrea & Phil | Test MaxRange in both azimuth and elevation. |
| 2025/03/21 14:13       | Andrea & Phil | Black mark shifted by 7-8 deg Indication that something has moved. Therefore, we added a new green mark on elevation = 80 deg. |
| 2025/04/01 11:10       | Andrea & Marco | We put the protection built by Philipp midway between the telescope and the satellite antenna. It has removed at the end of the day. |
| 2025/04/02 07:40       | Andrea & Marco | We put the protection built by Philipp just behind the telescope. |
| 2025/04/02 13:07       | Andrea & Marco | The protection has been removed. |
| 2025/04/02 15:00       | Andrea | According to the solar scans, I changed aziref to 181.2 |
| 2025/04/04 06:00       | Andrea | The change I've done on 2025/04/02 was wrong, it should have been done on the other direction! I changed aziref now to 182.4 |
| 2025/04/10 11:30       | Phil | The reflecting papers have been significantly reduced in size. We need to check now whether we still have the 15 minutes drift. |
| 2025/04/10 15:00       | Phil | Measurements and other stuff done at the telescope (dish + feed mainly). |
| 2025/04/22 08:46-09:49 | Phil | Measurements and other stuff done at the telescope (dish + feed mainly). |
| 2025/05/05 17:38       | Daniel/Marco  | Changed aziref from 182.4 to 183.4 |
| 2025/05/05 ~11:00      | Philipp/Maurizio  | Changed aziref from 182.4 to 186.3 (65.89° elevation) |