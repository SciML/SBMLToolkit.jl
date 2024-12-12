function is_event_assignment(k, model)
    for ev in last.(model.events)
        for as in ev.event_assignments
            if as.variable == k
                return true
            end
        end
    end
    return false
end
