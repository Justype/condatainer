-- Lmod modulefile.lua
-- Author: Ziyue Cheng
-- Lmod modulefile for {APP_NAME} {APP_VERSION}
-- Almost all of them are automated, you can modify whatis, ModulesHelp as needed.
-- You can also load other dependency modules by depends_on().
--     depends_on("openjdk/17.0.12")

local abs_path = myFileName()                   -- absolute path of module, /somewhere/modules/apps-modules/app/x.y.z
local app_full_name = myModuleFullName()        -- full name of module, app/x.y.z
local app_name, app_version = app_full_name:match("([^/]+)/([^/]+)")
-- local app_version = myModuleVersion()
local module_root = abs_path
for i = 1, 3 do
    module_root = module_root:match("(.+)/[^/]+$")
end    -- abs path of modules, /somewhere/modules
local app_root = pathJoin(module_root, "apps", app_full_name) -- abs path of target app, /somewhere/modules/apps/app/x.y.z
local current_datetime = os.date("%Y-%m-%d %H:%M:%S") -- Record current datetime

whatis("{WHATIS}")
help("{HELP}")

-- if (mode() == "load") then
--     io.stderr:write("[" .. current_datetime .. "] Loading module " .. app_full_name .. "\n")
-- -- elseif (mode() == "unload") then
-- end

if (isDir(pathJoin(app_root, "bin"))) then
    prepend_path("PATH", pathJoin(app_root, "bin"))
elseif (mode() == "load" and os.getenv("SLURM_JOB_ID") == nil) then
    io.stderr:write("WARNING: No bin directory found in " .. app_full_name .. "\n")
end

local env_var_name = string.upper(app_name:gsub("[-%./]", "_")) .. "_HOME"
setenv(env_var_name, app_root)

-- Handle conflicts with other versions
conflict(app_name)

