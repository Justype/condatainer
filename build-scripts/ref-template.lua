-- Lmod modulefile.lua
-- Author: Ziyue Cheng
-- Lmod modulefile for reference data
-- Almost all of them are automated, you can modify whatis, ModulesHelp as needed.

local abs_path = myFileName()                   -- absolute path of module, /somewhere/modules/ref-modules/assembly/data-type/version
local ref_full_name = myModuleFullName()        -- full name of module, assembly/data-type/version
local assembly, data_type, version = ref_full_name:match("([^/]+)/([^/]+)/([^/]+)")
local module_root = abs_path
for i = 1, 4 do
    module_root = module_root:match("(.+)/[^/]+$")
end    -- abs path of modules, /somewhere/modules
local ref_root = pathJoin(module_root, "ref", ref_full_name) -- abs path of target ref
local app_root = ref_root -- Make this script compatible with app modulefiles
local current_datetime = os.date("%Y-%m-%d %H:%M:%S") -- Record current datetime

whatis("Loads " .. ref_full_name)
help("Assembly: " .. assembly .. "\tData type: " .. data_type .. "\tVersion: " .. version .. "")

-- if (mode() == "load") then
--     io.stderr:write("[" .. current_datetime .. "] Loading module " .. ref_full_name .. "\n")
-- -- elseif (mode() == "unload") then
-- end

local env_var_name = string.upper(ref_full_name:gsub("[-%./]", "_")) .. "_HOME"
setenv(env_var_name, ref_root)

