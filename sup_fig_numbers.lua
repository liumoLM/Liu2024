-- Quarto numbers figures as "<fig-title> <n>", so with fig-title "Figure S"
-- the captions read "Figure S 1". Remove the space to give "Figure S1".
-- Runs after Quarto has written the caption labels (at: post-render).
local function is_gap(el)
  return el.t == "Space" or (el.t == "Str" and el.text == "\u{a0}")
end
function Inlines(inl)
  local i = 1
  while i <= #inl - 2 do
    if inl[i].t == "Str" and inl[i].text:match("S$") and is_gap(inl[i + 1])
        and inl[i + 2].t == "Str" and inl[i + 2].text:match("^%d") then
      inl[i].text = inl[i].text .. inl[i + 2].text
      inl:remove(i + 2)
      inl:remove(i + 1)
    else
      i = i + 1
    end
  end
  -- Labels may also arrive as a single Str with a no-break space.
  for _, el in ipairs(inl) do
    if el.t == "Str" then
      el.text = el.text:gsub("S\u{a0}(%d)", "S%1")
    end
  end
  return inl
end
