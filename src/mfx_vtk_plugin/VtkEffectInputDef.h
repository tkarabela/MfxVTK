#pragma once

#include <string>
#include <vector>
#include <PluginSupport/MfxEffect>

struct VtkEffectAttributeDef {
    MfxAttributeAttachment attachment;
    const char *name; // non-owning reference; preferably use static strings
    int componentCount;
    MfxAttributeType type;
    MfxAttributeSemantic semantic;
    bool mandatory;
};

struct VtkEffectInputDef {
    const char *name; // non-owning reference; preferably use static strings
    const char *label; // non-owning reference; preferably use static strings
    bool request_geometry;
    bool request_transform;
    bool is_output;
    std::vector<VtkEffectAttributeDef> requested_attributes;

    VtkEffectInputDef(const char *label, bool is_output=false);

    VtkEffectInputDef & Label(const char *label);
    VtkEffectInputDef & RequestAttribute(MfxAttributeAttachment attachment, const char* name, int componentCount, MfxAttributeType type, MfxAttributeSemantic semantic, bool mandatory);
    VtkEffectInputDef & RequestPointAttribute(const char* name, int componentCount, MfxAttributeType type, MfxAttributeSemantic semantic, bool mandatory);
    VtkEffectInputDef & RequestVertexAttribute(const char* name, int componentCount, MfxAttributeType type, MfxAttributeSemantic semantic, bool mandatory);
    VtkEffectInputDef & RequestFaceAttribute(const char* name, int componentCount, MfxAttributeType type, MfxAttributeSemantic semantic, bool mandatory);
    VtkEffectInputDef & RequestMeshAttribute(const char* name, int componentCount, MfxAttributeType type, MfxAttributeSemantic semantic, bool mandatory);
    VtkEffectInputDef & RequestGeometry(bool request);
    VtkEffectInputDef & RequestTransform(bool request);
};
