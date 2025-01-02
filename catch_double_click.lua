-- Back button is double-clicking.
-- Unbind it in Assignments, create Back macro that calls Back mouse cmd
-- Logitech G HUB > MANAGE PROFILES > Scripting

lastInstance = 0

function OnEvent(event, arg)
if event == "MOUSE_BUTTON_PRESSED" and arg == 4 then
    if (GetRunningTime() - lastInstance) >= 250 then
        PlayMacro("Back")
    end
    lastInstance = GetRunningTime();
end
end
