local keywords = {"Float", "Int", "Vector", "Geometry", "Bool", "Matrix", "Rotation", "Material", "Color", "Collection", "String", "Name", "Object"}

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
end

