function bodyFrame = getOsimFrame(OsimModel, bodyName)

    % get frame list
    frameList = OsimModel.getFrameList();

    % select frame based on the body name
    iter = frameList.begin();

    while 1   % run through all frames until get the same body name
        if iter.getName() == bodyName
            bodyFrame = iter.findBaseFrame();
            break;
        end

        iter.next();
    end
end