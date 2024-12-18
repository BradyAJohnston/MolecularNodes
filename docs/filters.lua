local keywords = {"Float", "Int", "Vector", "Geometry", "Bool", "Matrix", "Rotation", "Material", "Color", "Collection", "String", "Name", "Object", "Input"}
local attributes = {"sec_struct", "res_id", "chain_id", "entity_id", "res_name", "atom_name", "atom_id", "occupancy", "bfactor", "atomic_number", "charge", "vdw_radii", "mass", "is_backbone", "is_side_chain", "is_alpha_carbon", "Position", "Index", "Normal", "ID"}

local combined = {}
for _, v in ipairs(keywords) do
  table.insert(combined, v)
end
for _, v in ipairs(attributes) do
  table.insert(combined, v)
end

function Code(el)
  for _, keyword in ipairs(keywords) do
    local pattern = "(.+)::" .. keyword
    local name = el.text:match(pattern)
    if name then
      table.insert(el.classes, "custom-" .. keyword:lower())
      el.text = name
      return el
    end
  end
  for _, attribute in ipairs(combined) do
    if el.text == attribute then
      -- el.text = "[".. el.text .. "](attributes.qmd#" .. attribute .. ")"
      table.insert(el.classes, "custom-attribute")
      return el
    end
  end
end

