# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

prettyname(obj) = prettyname(typeof(obj))
function prettyname(T::Type)
  name = string(T)
  name = replace(name, r"{.*" => "")
  replace(name, r".+\." => "")
end
