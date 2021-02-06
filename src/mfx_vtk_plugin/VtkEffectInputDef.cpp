#include "VtkEffectInputDef.h"
#include <PluginSupport/MfxEffect>


VtkEffectInputDef::VtkEffectInputDef(const char *name, bool is_output)
        : name(name), request_geometry(true), request_transform(true), is_output(is_output), label(nullptr)
{}

VtkEffectInputDef & VtkEffectInputDef::Label(const char *label) {
    this->label = label;
    return *this;
}

VtkEffectInputDef & VtkEffectInputDef::RequestAttribute(MfxAttributeAttachment attachment, const char* name, int componentCount, MfxAttributeType type, MfxAttributeSemantic semantic, bool mandatory) {
    printf("MfxVTK - requesting attribute '%s'\n", name);
    requested_attributes.push_back(VtkEffectAttributeDef{
        .attachment = attachment,
        .name = name,
        .componentCount = componentCount,
        .type = type,
        .semantic = semantic,
        .mandatory = mandatory
    });
    return *this;
}

VtkEffectInputDef & VtkEffectInputDef::RequestPointAttribute(const char* name, int componentCount, MfxAttributeType type, MfxAttributeSemantic semantic, bool mandatory) {
    return RequestAttribute(MfxAttributeAttachment::Point, name, componentCount, type, semantic, mandatory);
}

VtkEffectInputDef & VtkEffectInputDef::RequestVertexAttribute(const char* name, int componentCount, MfxAttributeType type, MfxAttributeSemantic semantic, bool mandatory) {
    return RequestAttribute(MfxAttributeAttachment::Vertex, name, componentCount, type, semantic, mandatory);
}

VtkEffectInputDef & VtkEffectInputDef::RequestFaceAttribute(const char* name, int componentCount, MfxAttributeType type, MfxAttributeSemantic semantic, bool mandatory) {
    return RequestAttribute(MfxAttributeAttachment::Face, name, componentCount, type, semantic, mandatory);
}

VtkEffectInputDef & VtkEffectInputDef::RequestMeshAttribute(const char* name, int componentCount, MfxAttributeType type, MfxAttributeSemantic semantic, bool mandatory) {
    return RequestAttribute(MfxAttributeAttachment::Mesh, name, componentCount, type, semantic, mandatory);
}

VtkEffectInputDef & VtkEffectInputDef::RequestGeometry(bool request) {
    request_geometry = request;
    return *this;
}

VtkEffectInputDef & VtkEffectInputDef::RequestTransform(bool request) {
    request_transform = request;
    return *this;
}
