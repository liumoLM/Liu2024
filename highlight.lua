-- Map span classes pink / yellow / blue to background highlights.
local latex = { pink = "pink!45", yellow = "yellow!55", blue = "cyan!30" }
local html  = { pink = "#ffc0cb", yellow = "#ffff99", blue = "#b3ecff" }
function Span(el)
  for cls, col in pairs(latex) do
    if el.classes:includes(cls) then
      if FORMAT:match("latex") then
        local out = pandoc.Inlines({ pandoc.RawInline("latex", "\\highLight[" .. col .. "]{") })
        out:extend(el.content)
        out:insert(pandoc.RawInline("latex", "}"))
        return out
      elseif FORMAT:match("html") then
        el.attributes["style"] = "background-color:" .. html[cls]
        return el
      end
    end
  end
end
