%% Old Datetime Extraction
% Pre-background measurements
    % Pre-Background
        VCAC_11Jul2017_PreBack_DT     = datetime(['11-Jul-2017 11:09:17';
                                        '11-Jul-2017 11:44:59';
                                        '11-Jul-2017 12:27:15';
                                        '11-Jul-2017 12:54:44';
                                        '11-Jul-2017 11:13:47';
                                        '11-Jul-2017 11:51:39';
                                        '11-Jul-2017 12:31:58';
                                        '11-Jul-2017 12:59:40']);
       VCAC_11Jul2017_PreBack_DT      = reshape(VCAC_11Jul2017_PreBack_DT, [nchams 2]);

    % Chamber ON
        VCAC_11Jul2017_ChamON_DT      = datetime(['11-Jul-2017 11:13:51';
                                        '11-Jul-2017 11:51:43';
                                        '11-Jul-2017 12:32:01';
                                        '11-Jul-2017 13:00:29';
                                        '11-Jul-2017 11:22:24';
                                        '11-Jul-2017 11:58:55';
                                        '11-Jul-2017 12:36:56';
                                        '11-Jul-2017 13:10:55']);
        VCAC_11Jul2017_ChamON_DT      = reshape(VCAC_11Jul2017_ChamON_DT, [nchams 2]);

    % Post-Background
        VCAC_11Jul2017_PostBack_DT    = datetime(['11-Jul-2017 11:22:28';
                                        '11-Jul-2017 11:58:58';
                                        '11-Jul-2017 12:37:29';
                                        '11-Jul-2017 13:10:20';
                                        '11-Jul-2017 11:28:23';
                                        '11-Jul-2017 12:04:59';
                                        '11-Jul-2017 12:42:30';
                                        '11-Jul-2017 13:15:19']);
        VCAC_11Jul2017_PostBack_DT    = reshape(VCAC_11Jul2017_PostBack_DT, [nchams 2]);

    % Date Serials
        % Pre-Background
            VCAC_11Jul2017_PreBack_DS = datenum(VCAC_11Jul2017_PreBack_DT);
        % Chamber ON
            VCAC_11Jul2017_ChamON_DS  = datenum(VCAC_11Jul2017_ChamON_DT);
        % Post-Background
            VCAC_11Jul2017_PostBack_DS= datenum(VCAC_11Jul2017_PostBack_DT);
