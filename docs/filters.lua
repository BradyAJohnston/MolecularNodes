local keywords = {"Float", "Int", "Vector", "Geometry", "Bool", "Matrix", "Rotation", "Material", "Color", "Collection", "String"}

function Code(el)
  for _, keyword in ipairs(keywords) do
    if el.text == keyword then
      table.insert(el.classes, "custom-" .. keyword:lower())
      return el
    end
  end
end

-- function Str(el)
--   for _, keyword in ipairs(keywords) do
--     if el.text == keyword then
--       return pandoc.Span(el.text, pandoc.Attr("", {"custom-" .. keyword:lower()}, {}))
--     end
--   end
-- end